"""The public filter of the docs site and the GitHub wiki mirror.

The experiments wiki is the author's working record and names the machines a run used, the
local libraries and the private repository. None of that belongs on a public page. Every
Markdown page that leaves for the docs site or the GitHub wiki passes through filter_markdown()
and is then audited: a page that still carries a private token fails the build.

    python3 site/public.py --filter <src.md> <dst.md>     filter one page (exit 1 on a residue)
    python3 site/public.py --audit <file> ...             report residues, exit 1 if any
"""
import re
import sys

# rows of the status table at the head of an experiment page that are not public
DROP_ROWS = ("When", "Where")
# the private repository on the forge: links to its files become the file path in code
FORGE_SRC = re.compile(r"\[([^\]]+)\]\(https://whitebox\.tailf2ef50\.ts\.net:3000/pancho/[^/)]+/src/branch/[^/)]+/([^)\s]+)\)")
FORGE_ANY = re.compile(r"\[([^\]]+)\]\(https://whitebox\.tailf2ef50\.ts\.net:3000/[^)\s]+\)")
# what must not appear on a public page (machine names, the tailnet, local paths, libraries)
PRIVATE = [
    (re.compile(r"tailf2ef50|tailscale", re.I), "the tailnet address"),
    (re.compile(r"\bwhitebox\b", re.I), "a machine name (whitebox)"),
    (re.compile(r"\bStudio\b"), "a machine name (Studio)"),
    (re.compile(r"\b(wbx|pro2018|panchos-ia|statisticals)\b", re.I), "a machine name"),
    (re.compile(r"(^|[\s(`'\"])~/"), "a home-directory path"),
    (re.compile(r"/Users/|/home/"), "a local path"),
    (re.compile(r"\brlib-"), "a local R library"),
    (re.compile(r"\bssh\b"), "a shell command on a private host"),
]


def filter_markdown(text):
    out = []
    for line in text.split("\n"):
        m = re.match(r"^\|\s*([^|]+?)\s*\|", line)
        if m and m.group(1) in DROP_ROWS:
            continue
        line = FORGE_SRC.sub(lambda m: "`%s`" % m.group(2), line)
        line = FORGE_ANY.sub(lambda m: m.group(1), line)
        out.append(line)
    return "\n".join(out)


def audit(text, name="<text>"):
    found = []
    for i, line in enumerate(text.split("\n"), 1):
        for rx, what in PRIVATE:
            if rx.search(line):
                found.append("%s:%d: %s: %s" % (name, i, what, line.strip()[:120]))
                break
    return found


def main(argv):
    if argv[:1] == ["--filter"] and len(argv) == 3:
        text = filter_markdown(open(argv[1], encoding="utf-8").read())
        open(argv[2], "w", encoding="utf-8").write(text)
        bad = audit(text, argv[2])
        print("\n".join(bad))
        return 1 if bad else 0
    if argv[:1] == ["--audit"]:
        bad = []
        for f in argv[1:]:
            bad += audit(open(f, encoding="utf-8").read(), f)
        print("\n".join(bad) if bad else "public audit: clean (%d files)" % (len(argv) - 1))
        return 1 if bad else 0
    print(__doc__)
    return 2


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
