#!/usr/bin/env bash
#
# publish.sh -- share the SHUTTLE project with the team via the COS bucket.
#
# Every run produces TWO AES-256 encrypted zip archives (extraction password
# below) and uploads them into the bucket's SHUTTLE-Code/ directory:
#
#   SHUTTLE-<TS>.zip       the whole SHUTTLE/ tree (minus .git, .github,
#                          publish.sh, and the regenerable dist/ build dir)
#   SHUTTLE-NGCC-<TS>.zip  the NGCC submission products that pack_ngcc.sh
#                          generates (Implementations/ + Test_Vectors/)
#
# where <TS> is a yyyymmddHHMM timestamp (24-hour), e.g. 202606251439.
#
# After the two new archives are uploaded, every older top-level package in
# SHUTTLE-Code/ is moved into SHUTTLE-Code/Archive/, and SHUTTLE-Code/index.html
# is rebuilt: a collapsible listing of everything under SHUTTLE-Code/ (the
# current archives at the top, the Archive/ folder collapsed by default).
#
# The decryption password is shared with the team out of band; it is NOT shown
# in index.html. This publish.sh is excluded from the SHUTTLE-<TS>.zip snapshot.
# NOTE: 7z takes the password as a command-line arg (-p...), briefly visible to
# other local users via 'ps'; this script is meant for single-user hosts.
#
# Usage:
#   ./publish.sh                 build + upload both archives, rotate, reindex
#   ./publish.sh --no-upload     build both archives locally only (no COS calls)
#   ./publish.sh --out DIR       with --no-upload, where to drop the archives
#   ./publish.sh --no-sanity     skip pack_ngcc.sh's standalone-build/KAT check
#   ./publish.sh --refresh-index only rebuild + upload SHUTTLE-Code/index.html
#   ./publish.sh -h | --help
#
# Environment overrides:
#   SHUTTLE_ZIP_PASSWORD  extraction password (default: the embedded one)
#   COSCMD_BIN            path to coscmd
#   SEVENZIP_BIN          path to 7z / 7za / 7zz
#   PROXY_URL             public URL prefix (default https://share.ji-peng.com)
#   COS_HOST             COS host for server-side move (default derived from
#                         ~/.cos.conf as <bucket>.cos.<region>.myqcloud.com)

set -euo pipefail

# --- configuration ----------------------------------------------------------

BUCKET_DIR="SHUTTLE-Code"
ARCHIVE_DIR="${BUCKET_DIR}/Archive"
INDEX_KEY="${BUCKET_DIR}/index.html"
PROXY_URL="${PROXY_URL:-https://share.ji-peng.com}"
COS_CONF="${COS_CONF:-$HOME/.cos.conf}"
ZIP_PASSWORD="${SHUTTLE_ZIP_PASSWORD:-hycFkBh+okdvxYX2c5vbOwJGR7fTg/DZ}"

SCRIPT_PATH="$(readlink -f "${BASH_SOURCE[0]}")"
SCRIPT_DIR="$(dirname "$SCRIPT_PATH")"   # = the SHUTTLE/ project root

# --- small helpers ----------------------------------------------------------

die() {
  echo "Error: $*" >&2
  exit 1
}

resolve_coscmd() {
  if [[ -n "${COSCMD_BIN:-}" ]]; then
    [[ -x "$COSCMD_BIN" ]] || die "COSCMD_BIN is set but not executable: $COSCMD_BIN"
    printf '%s\n' "$COSCMD_BIN"
    return 0
  fi
  command -v coscmd >/dev/null 2>&1 || die "coscmd not found. Set COSCMD_BIN or install coscmd."
  command -v coscmd
}

resolve_7z() {
  if [[ -n "${SEVENZIP_BIN:-}" ]]; then
    [[ -x "$SEVENZIP_BIN" ]] || die "SEVENZIP_BIN is set but not executable: $SEVENZIP_BIN"
    printf '%s\n' "$SEVENZIP_BIN"
    return 0
  fi
  local cand
  for cand in 7z 7zz 7za; do
    if command -v "$cand" >/dev/null 2>&1; then
      command -v "$cand"
      return 0
    fi
  done
  die "7z not found (need 7z, 7zz, or 7za). Set SEVENZIP_BIN or install p7zip."
}

# Read a non-secret key from the INI-style ~/.cos.conf (bucket / region).
cos_conf_get() {
  local key="$1"
  [[ -f "$COS_CONF" ]] || return 1
  awk -F= -v k="$key" '
    { line=$0; sub(/^[ \t]+/, "", line) }
    line ~ "^"k"[ \t]*=" {
      val=$2; gsub(/[ \t\r]/, "", val); print val; exit
    }
  ' "$COS_CONF"
}

# Host form coscmd move/copy expects: <bucket>.cos.<region>.myqcloud.com
resolve_cos_host() {
  if [[ -n "${COS_HOST:-}" ]]; then
    printf '%s\n' "$COS_HOST"
    return 0
  fi
  local bucket region
  bucket="$(cos_conf_get bucket || true)"
  region="$(cos_conf_get region || true)"
  [[ -n "$bucket" && -n "$region" ]] || \
    die "could not derive COS host from $COS_CONF (need bucket + region); set COS_HOST."
  printf '%s.cos.%s.myqcloud.com\n' "$bucket" "$region"
}

# Create a tzip archive with AES-256 encryption. Args after the archive path
# are passed to 7z verbatim (inputs and -xr!/-x! excludes). NOTE: no '--'
# separator -- it would stop switch parsing and turn '-xr!.git' into a literal
# filename. 7z exit codes: 0 ok, 1 non-fatal warning, >=2 fatal.
make_zip() {
  local archive="$1"; shift
  local rc=0
  "$SEVENZIP" a -tzip -mx=5 -mem=AES256 -p"$ZIP_PASSWORD" -bso0 -bsp0 \
    "$archive" "$@" || rc=$?
  (( rc <= 1 )) || die "7z failed (exit $rc) building $archive"
  [[ -s "$archive" ]] || die "7z produced an empty archive: $archive"
}

# --- build the two archives -------------------------------------------------
# $1 = output directory for the archives. Sets SHUTTLE_ZIP / NGCC_ZIP / NGCC_BUILD.

build_archives() {
  local outdir="$1"
  mkdir -p "$outdir"

  SHUTTLE_ZIP="${outdir}/SHUTTLE-${TS}.zip"
  NGCC_ZIP="${outdir}/SHUTTLE-NGCC-${TS}.zip"

  # 1) Whole SHUTTLE/ tree. Build from inside the project so paths are
  #    relative to the project root. Exclude VCS dirs, this script, and the
  #    regenerable pack_ngcc.sh output dir (its contents are the NGCC zip).
  echo "Packing ${SHUTTLE_ZIP##*/} (SHUTTLE/ source snapshot) ..."
  rm -f "$SHUTTLE_ZIP"
  (
    cd "$SCRIPT_DIR"
    make_zip "$SHUTTLE_ZIP" . \
      '-xr!.git' '-xr!.github' '-xr!publish.sh' '-xr!dist'
  )

  # 2) NGCC submission products. Generate them into a throwaway dir OUTSIDE
  #    SHUTTLE/ (so they never pollute the source snapshot above), then zip
  #    the Implementations/ + Test_Vectors/ trees at the archive root.
  echo "Generating NGCC submission products (pack_ngcc.sh) ..."
  NGCC_BUILD="$(mktemp -d)"
  TMP_DIRS+=("$NGCC_BUILD")
  local pack_args=(--out "$NGCC_BUILD")
  (( NO_SANITY )) && pack_args+=(--no-sanity)
  bash "${SCRIPT_DIR}/pack_ngcc.sh" "${pack_args[@]}"
  [[ -d "${NGCC_BUILD}/Implementations" && -d "${NGCC_BUILD}/Test_Vectors" ]] || \
    die "pack_ngcc.sh did not produce Implementations/ + Test_Vectors/ in $NGCC_BUILD"

  echo "Packing ${NGCC_ZIP##*/} (NGCC submission products) ..."
  rm -f "$NGCC_ZIP"
  (
    cd "$NGCC_BUILD"
    local items=(Implementations Test_Vectors)
    [[ -e README ]] && items+=(README)
    make_zip "$NGCC_ZIP" "${items[@]}"
  )
}

# --- COS interaction --------------------------------------------------------

upload_archives() {
  echo "Uploading archives to ${BUCKET_DIR}/ ..."
  "$COSCMD" upload "$SHUTTLE_ZIP" "${BUCKET_DIR}/$(basename "$SHUTTLE_ZIP")"
  "$COSCMD" upload "$NGCC_ZIP"   "${BUCKET_DIR}/$(basename "$NGCC_ZIP")"
  echo "Uploaded:"
  echo "  ${PROXY_URL}/${BUCKET_DIR}/$(basename "$SHUTTLE_ZIP")"
  echo "  ${PROXY_URL}/${BUCKET_DIR}/$(basename "$NGCC_ZIP")"
}

# Move every older top-level package in SHUTTLE-Code/ into SHUTTLE-Code/Archive/.
# "Top-level" = directly under SHUTTLE-Code/ (no further '/'); the two archives
# just uploaded this run are kept in place.
rotate_archives() {
  local keep_a keep_b host listing rel name
  keep_a="$(basename "$SHUTTLE_ZIP")"
  keep_b="$(basename "$NGCC_ZIP")"
  host="$(resolve_cos_host)"

  listing="$("$COSCMD" list "${BUCKET_DIR}/" -r 2>/dev/null || true)"
  [[ -n "$listing" ]] || return 0

  local moved=0
  while IFS= read -r key; do
    [[ -n "$key" ]] || continue
    rel="${key#"${BUCKET_DIR}/"}"          # strip the SHUTTLE-Code/ prefix
    [[ "$rel" == "$key" ]] && continue     # not under SHUTTLE-Code/, skip
    [[ "$rel" == */* ]] && continue        # in a subdir (e.g. Archive/), skip
    [[ "$rel" == *.zip ]] || continue      # only archive packages
    name="$rel"
    [[ "$name" == "$keep_a" || "$name" == "$keep_b" ]] && continue  # keep new
    echo "Archiving older package: $name"
    "$COSCMD" move "${host}/${BUCKET_DIR}/${name}" "${ARCHIVE_DIR}/${name}"
    moved=$((moved + 1))
  done < <(printf '%s\n' "$listing" | awk 'NF>=5 && $2 ~ /^[0-9]+$/ {print $1}')

  echo "Archived ${moved} older package(s)."
}

# Rebuild SHUTTLE-Code/index.html: a collapsible tree of everything under
# SHUTTLE-Code/. Top-level files listed directly; each subdirectory (Archive/)
# rendered as a <details> collapsed by default.
update_index() {
  command -v python3 >/dev/null 2>&1 || {
    echo "Warning: python3 not found; skipping index refresh." >&2
    return 0
  }
  echo "Rebuilding ${INDEX_KEY} ..."

  local tmp_list tmp_index
  tmp_list="$(mktemp)"
  tmp_index="$(mktemp --suffix=.html)"
  TMP_FILES+=("$tmp_list" "$tmp_index")

  if ! "$COSCMD" list "${BUCKET_DIR}/" -r >"$tmp_list" 2>/dev/null; then
    echo "Warning: failed to list ${BUCKET_DIR}/; index.html not refreshed." >&2
    return 0
  fi

  if ! PROXY_URL="$PROXY_URL" BUCKET_DIR="$BUCKET_DIR" \
       python3 - "$tmp_list" "$tmp_index" <<'PY'
import sys, os, re, html, datetime
from collections import defaultdict

list_path, out_path = sys.argv[1], sys.argv[2]
proxy = os.environ.get("PROXY_URL", "").rstrip("/")
base = os.environ.get("BUCKET_DIR", "").strip("/")
prefix = base + "/"
date_re = re.compile(r"^\d{4}-\d{2}-\d{2}$")

# key -> (size_bytes, "YYYY-MM-DD HH:MM:SS")
top_files = []                  # (name, key, size, when)
sub = defaultdict(list)         # subdir -> [(name, key, size, when)]

def human(n):
    f = float(n)
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if f < 1024 or unit == "TB":
            return (f"{int(f)} {unit}" if unit == "B" else f"{f:.1f} {unit}")
        f /= 1024

with open(list_path, encoding="utf-8", errors="replace") as f:
    for line in f:
        fields = line.split()
        if len(fields) < 5 or not fields[1].isdigit() or not date_re.match(fields[-2]):
            continue
        key, size = fields[0], int(fields[1])
        when = fields[-2] + " " + fields[-1]
        if key.endswith("/"):
            continue
        if not key.startswith(prefix):
            continue
        rel = key[len(prefix):]
        if rel in ("", "index.html"):
            continue
        if "/" in rel:
            folder, _, name = rel.partition("/")
            sub[folder].append((name, key, size, when))
        else:
            top_files.append((rel, key, size, when))

now = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
total = len(top_files) + sum(len(v) for v in sub.values())

def esc(s):
    return html.escape(str(s), quote=True)

def row(name, key, size, when):
    url = f"{proxy}/{key}"
    return (f"<li><a href=\"{esc(url)}\" target=\"_blank\" "
            f"rel=\"noopener noreferrer\">{esc(name)}</a>"
            f"<span class=\"meta\">{esc(human(size))} &middot; {esc(when)}</span></li>")

parts = [
    "<!doctype html>",
    "<html lang=\"zh-CN\">",
    "<head>",
    "<meta charset=\"utf-8\">",
    "<meta name=\"viewport\" content=\"width=device-width,initial-scale=1\">",
    "<title>SHUTTLE-Code</title>",
    "<style>",
    "body{font-family:-apple-system,BlinkMacSystemFont,\"Segoe UI\",Helvetica,Arial,sans-serif;"
    "max-width:920px;margin:2rem auto;padding:0 1rem;color:#222;line-height:1.6;}",
    "h1{border-bottom:2px solid #eee;padding-bottom:.4rem;margin-bottom:.4rem;}",
    "ul{list-style:none;padding-left:1rem;margin:.3rem 0;}",
    "li{padding:.15rem 0;}",
    "a{color:#0366d6;text-decoration:none;word-break:break-all;}",
    "a:hover{text-decoration:underline;}",
    ".meta{color:#888;font-size:.85em;margin-left:.6rem;}",
    ".note{color:#888;font-size:.9em;margin-bottom:1rem;}",
    "details{margin:.3rem 0;}",
    "summary{cursor:pointer;font-weight:600;color:#333;}",
    ".empty{color:#888;}",
    "</style>",
    "</head>",
    "<body>",
    "<h1>SHUTTLE-Code</h1>",
    f"<p class=\"note\">Updated at {esc(now)} &middot; {total} file(s) &middot; "
    "archives are AES-256 encrypted (password shared separately).</p>",
]

if not top_files and not sub:
    parts.append("<p class=\"empty\">No files found.</p>")
else:
    parts.append("<ul>")
    for name, key, size, when in sorted(top_files, key=lambda x: x[0].lower(), reverse=True):
        parts.append("  " + row(name, key, size, when))
    parts.append("</ul>")
    for folder in sorted(sub.keys(), key=str.lower):
        items = sorted(sub[folder], key=lambda x: x[0].lower(), reverse=True)
        parts.append(f"<details><summary>{esc(folder)}/ ({len(items)})</summary>")
        parts.append("<ul>")
        for name, key, size, when in items:
            parts.append("  " + row(name, key, size, when))
        parts.append("</ul></details>")

parts.append("</body></html>")
with open(out_path, "w", encoding="utf-8") as f:
    f.write("\n".join(parts) + "\n")
PY
  then
    echo "Warning: failed to build index.html." >&2
    return 0
  fi

  [[ -s "$tmp_index" ]] || { echo "Warning: empty index.html; skipping." >&2; return 0; }

  "$COSCMD" upload \
    -H "Cache-Control: no-cache, no-store, must-revalidate" \
    "$tmp_index" "$INDEX_KEY" >/dev/null
  echo "Index refreshed: ${PROXY_URL}/${INDEX_KEY}"
}

# --- main -------------------------------------------------------------------

main() {
  local no_upload=0 refresh_only=0 out_dir=""
  NO_SANITY=0
  TMP_DIRS=()
  TMP_FILES=()

  while (( $# > 0 )); do
    case "$1" in
      --no-upload)     no_upload=1; shift ;;
      --no-sanity)     NO_SANITY=1; shift ;;
      --refresh-index) refresh_only=1; shift ;;
      --out)           shift; (( $# > 0 )) || die "--out requires a directory."; out_dir="$1"; shift ;;
      --out=*)         out_dir="${1#*=}"; shift ;;
      -h|--help)       awk 'NR==1{next} /^#/{sub(/^# ?/,"");print;next} {exit}' "$SCRIPT_PATH"; return 0 ;;
      *)               die "unknown argument: $1 (try --help)" ;;
    esac
  done

  SEVENZIP="$(resolve_7z)"
  [[ -n "$ZIP_PASSWORD" ]] || \
    die "archive password is empty; refusing to create an unencrypted zip (set SHUTTLE_ZIP_PASSWORD)."
  TS="$(date +%Y%m%d%H%M)"

  # Clean up temp dirs/files on exit (built archives in a temp dir are removed;
  # with --no-upload we copy them to out_dir first, see below).
  cleanup() {
    local d f
    for d in "${TMP_DIRS[@]:-}"; do [[ -n "$d" && -d "$d" ]] && rm -rf "$d"; done
    for f in "${TMP_FILES[@]:-}"; do [[ -n "$f" && -f "$f" ]] && rm -f "$f"; done
    return 0   # EXIT-trap last status leaks to $?; force success on clean exit.
  }
  trap cleanup EXIT

  if (( refresh_only )); then
    COSCMD="$(resolve_coscmd)"
    update_index
    return 0
  fi

  if (( no_upload )); then
    # Build locally only; leave the archives where the user can inspect them.
    local dest="${out_dir:-$PWD}"
    mkdir -p "$dest"
    local work; work="$(mktemp -d)"; TMP_DIRS+=("$work")
    build_archives "$work"
    cp -f "$SHUTTLE_ZIP" "$NGCC_ZIP" "$dest"/
    echo "Built (no upload):"
    echo "  ${dest%/}/$(basename "$SHUTTLE_ZIP")"
    echo "  ${dest%/}/$(basename "$NGCC_ZIP")"
    return 0
  fi

  COSCMD="$(resolve_coscmd)"
  local work; work="$(mktemp -d)"; TMP_DIRS+=("$work")
  build_archives "$work"
  upload_archives
  rotate_archives
  update_index
  echo "Done."
}

main "$@"
