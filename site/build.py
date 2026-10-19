#!/usr/bin/env python3
"""Build the emphasis site: Markdown sources and the experiments wiki -> the root pages.

    python3 site/build.py --wiki <path to the wiki checkout> [--out <repository root>]

Inputs
    site/template.html          the one page template
    site/style.css              the one stylesheet (copied to assets/site.css)
    site/src/index.md           -> index.html
    site/src/usage.md           -> usage.html
    site/reference/*.html       Rd fragments (site/reference.R, run on whitebox) -> reference.html
                                ordered by site/reference/index.json when it exists
    <wiki>/What-we-know.md      -> established.html
    <wiki>/Estimability.md      -> estimability.html
    <wiki>/Real-phylogenies.md  -> real-phylogenies.html
    <wiki>/Open-issues.md       -> open-issues.html
    <wiki>/Glossary.md          -> glossary.html
    <wiki>/Home.md              -> experiments/index.html (the table of experiments, grouped as Home groups it)
    <wiki>/E<n>.md, S.md        -> experiments/<name>.html
    <wiki>/images/*.png         -> assets/wiki/

Markdown is converted by pandoc (-f gfm+tex_math_dollars -t html5 --mathjax); the maths is
rendered in the browser by the local MathJax (assets/mathjax/tex-svg.js). Wiki links
([E32](E32), [Home](Home), [Glossary](Glossary#the-models)) are rewritten to the site's pages;
links to the forge and to GitHub are left as they are. Images get a <figure> with a numbered
caption taken from the italic line under them. Pages with more than three h2/h3 headings get a
table of contents and numbered sections. Every internal link and image is checked at the end;
a broken one makes the build exit with status 1.
"""

import argparse
import html
import json
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

SITE = Path(__file__).resolve().parent
sys.path.insert(0, str(SITE))
import public   # the public filter and audit (site/public.py)
ROOT_DEFAULT = SITE.parent

GITHUB = "https://github.com/franciscorichter/emphasis"

NAV = [
    ("Overview", "index.html"),
    ("Usage", "usage.html"),
    ("Reference", "reference.html"),
    ("What is established", "established.html"),
    ("Estimability", "estimability.html"),
    ("Real phylogenies", "real-phylogenies.html"),
    ("Experiments", "experiments/index.html"),
    ("Glossary", "glossary.html"),
    ("Simulator", "simulator.html"),
    ("D", "D.html"),
    ("GitHub", GITHUB),
]

# Wiki page name -> site path (relative to the repository root).
WIKI_PAGES = {
    "Home": "experiments/index.html",
    "What-we-know": "established.html",
    "Estimability": "estimability.html",
    "Real-phylogenies": "real-phylogenies.html",
    "Open-issues": "open-issues.html",
    "Glossary": "glossary.html",
}

EXPERIMENT_RE = re.compile(r"^(E\d+|S)$")

# Sections of the site that are kept as they are and only linked.
KEPT = ["D.html", "simulator.html"]

# Anchors the kept pages link to, kept stable whatever the heading's own id becomes:
# page -> {anchor: regex matched against the heading id it should land on}.
STABLE_ANCHORS = {"index.html": {"covariates": r"covariates$"}}


# ----------------------------------------------------------------------------- helpers

def die(msg):
    sys.stderr.write("build.py: " + msg + "\n")
    sys.exit(1)


def read(p):
    return Path(p).read_text(encoding="utf-8")


def write(p, text):
    p = Path(p)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(text, encoding="utf-8")


def pandoc(markdown):
    """Markdown (GitHub flavour with $ maths) -> an HTML5 body fragment."""
    cmd = [
        "pandoc",
        "-f", "gfm+tex_math_dollars",
        "-t", "html5",
        "--mathjax",
        "--wrap=none",
    ]
    r = subprocess.run(cmd, input=markdown, capture_output=True, text=True, encoding="utf-8")
    if r.returncode != 0:
        die("pandoc failed: " + r.stderr)
    return r.stdout


def strip_tags(s):
    return html.unescape(re.sub(r"<[^>]+>", "", s)).strip()


def rel(target, page):
    """Path of `target` relative to the directory of `page` (both relative to the root)."""
    page_dir = os.path.dirname(page)
    return os.path.relpath(target, page_dir) if page_dir else target


def root_prefix(page):
    depth = page.count("/")
    return "../" * depth


# ----------------------------------------------------------------------------- transforms

def wiki_target(name):
    """The site path of a wiki page name, or None if it is not one."""
    if name in WIKI_PAGES:
        return WIKI_PAGES[name]
    if EXPERIMENT_RE.match(name):
        return "experiments/" + name + ".html"
    return None


def rewrite_links(body, page):
    """[E32](E32) -> E32.html etc., relative to `page`; images -> assets/wiki/."""

    def href(m):
        url = m.group(1)
        if re.match(r"^[a-z][a-z0-9+.-]*:", url) or url.startswith(("/", "#", "../")):
            return m.group(0)
        name, _, frag = url.partition("#")
        target = wiki_target(name)
        if target is None:
            return m.group(0)
        new = rel(target, page) + ("#" + frag if frag else "")
        return 'href="' + new + '"'

    def src(m):
        url = m.group(1)
        if url.startswith("images/"):
            return 'src="' + rel("assets/wiki/" + url[len("images/"):], page) + '"'
        return m.group(0)

    body = re.sub(r'href="([^"]*)"', href, body)
    body = re.sub(r'src="([^"]*)"', src, body)
    return body


FIG_RE = re.compile(
    r'<p><img src="([^"]+)" alt="([^"]*)"\s*/?></p>\s*(?:<p><em>(.*?)</em></p>)?',
    re.S,
)


def figures(body):
    """A lone image, with the italic paragraph under it as its caption, -> a numbered <figure>."""
    n = [0]

    def repl(m):
        n[0] += 1
        src, alt, cap = m.group(1), m.group(2), m.group(3)
        caption = cap.strip() if cap else html.escape(html.unescape(alt))
        return (
            '<figure id="figure-%d">\n<img src="%s" alt="%s" loading="lazy">\n'
            '<figcaption><strong>Figure %d.</strong> %s</figcaption>\n</figure>'
            % (n[0], src, alt, n[0], caption)
        )

    return FIG_RE.sub(repl, body)


def tables(body, meta_first=False):
    """Drop pandoc's colgroups, wrap each table for narrow screens, mark the first table as metadata."""
    body = re.sub(r"<colgroup>.*?</colgroup>\s*", "", body, flags=re.S)
    first = [True]

    def repl(m):
        cls = ""
        if meta_first and first[0]:
            cls = ' class="meta"'
        first[0] = False
        return '<div class="table-wrap">\n<table%s>' % cls

    body = re.sub(r"<table>", repl, body)
    body = body.replace("</table>", "</table>\n</div>")
    return body


def take_h1(body):
    m = re.search(r"<h1[^>]*>(.*?)</h1>\s*", body, flags=re.S)
    if not m:
        return None, body
    return m.group(1).strip(), body[: m.start()] + body[m.end():]


HEAD_RE = re.compile(r'<h([23]) id="([^"]+)"[^>]*>(.*?)</h\1>', re.S)


def toc_and_numbering(body, numbered_allowed=True, run_in=False):
    """A table of contents from h2/h3 when there are more than three of them.

    Sections are numbered (CSS counters) unless the h2 headings already carry numbers or the
    page declines numbering. With run_in (the reference page, where the h3s are function
    names) each h2 becomes one line: the group, then its entries separated by commas.
    Returns (toc html or '', numbered flag)."""
    heads = HEAD_RE.findall(body)
    if len(heads) <= 3:
        return "", False
    if run_in:
        groups = []
        for lvl, hid, text in heads:
            text = re.sub(r'<span class="ref-title">.*?</span>', "", text, flags=re.S)
            text = re.sub(r"\s+", " ", text).strip()
            link = '<a href="#%s">%s</a>' % (hid, text)
            if lvl == "2" or not groups:
                groups.append([link])
            else:
                groups[-1].append(link)
        lines = ["<p>%s: %s.</p>" % (g[0], ", ".join(g[1:])) if len(g) > 1 else "<p>%s</p>" % g[0]
                 for g in groups]
        toc = ('<nav class="toc toc-index" aria-label="Contents">\n<p class="toc-title">Contents</p>\n%s\n</nav>'
               % "\n".join(lines))
        return toc, False
    numbered = numbered_allowed and not any(
        re.match(r"^\s*\d+[.)]", strip_tags(t)) for lvl, _, t in heads if lvl == "2")
    tag = "ol" if numbered else "ul"
    out = []
    depth = 0  # 0: no item open, 1: an h2 item open, 2: its h3 list open
    for lvl, hid, text in heads:
        text = re.sub(r'<span class="ref-title">.*?</span>', "", text, flags=re.S)  # reference: topic only
        text = re.sub(r"\s+", " ", text).strip()
        link = '<li><a href="#%s">%s</a>' % (hid, text)
        if lvl == "2":
            if depth == 2:
                out.append("</%s>" % tag)
            if depth >= 1:
                out.append("</li>")
            out.append(link)
            depth = 1
        else:
            if depth == 0:
                out.append("<li>")
                depth = 1
            if depth == 1:
                out.append("<%s>" % tag)
                depth = 2
            out.append(link + "</li>")
    if depth == 2:
        out.append("</%s>" % tag)
    if depth >= 1:
        out.append("</li>")
    toc = (
        '<nav class="toc" aria-label="Contents">\n<p class="toc-title">Contents</p>\n<%s>\n%s\n</%s>\n</nav>'
        % (tag, "\n".join(out), tag)
    )
    return toc, numbered


# ----------------------------------------------------------------------------- page assembly

class Site:
    def __init__(self, root, wiki, template, version):
        self.root = Path(root)
        self.wiki = Path(wiki) if wiki else None
        self.template = template
        self.version = version
        self.pages = []  # paths written, relative to root
        self.notes = []

    def nav(self, page):
        pre = root_prefix(page)
        items = []
        for label, target in NAV:
            if target.startswith("http"):
                items.append('<li><a href="%s" rel="external">%s</a></li>' % (target, label))
                continue
            active = page == target or (target == "experiments/index.html" and page.startswith("experiments/"))
            cur = ' aria-current="page"' if active else ""
            items.append('<li><a href="%s%s"%s>%s</a></li>' % (pre, target, cur, label))
        return "\n".join(items)

    def render(self, page, title, body, description="", meta_first=False, numbered_allowed=True, run_in_toc=False):
        """Assemble one page from a body fragment and write it."""
        body = rewrite_links(body, page)
        body = figures(body)
        body = tables(body, meta_first=meta_first)
        h1, body = take_h1(body)
        if h1 is None:
            h1 = html.escape(title)
        for anchor, pattern in STABLE_ANCHORS.get(page, {}).items():
            for m in re.finditer(r'<h[23] id="([^"]*)"', body):
                if re.search(pattern, m.group(1)):
                    body = body[: m.start()] + '<span id="%s"></span>' % anchor + body[m.start():]
                    break
            else:
                self.notes.append("%s: no heading for the stable anchor #%s" % (page, anchor))
        toc, numbered = toc_and_numbering(body, numbered_allowed, run_in=run_in_toc)
        pre = root_prefix(page)
        out = self.template
        subs = {
            "{{root}}": pre,
            "{{title}}": html.escape(strip_tags(h1)),
            "{{description}}": html.escape(description),
            "{{nav}}": self.nav(page),
            "{{h1}}": h1,
            "{{toc}}": toc,
            "{{body}}": body,
            "{{article_class}}": "numbered" if numbered else "",
            "{{version}}": self.version,
        }
        for k, v in subs.items():
            out = out.replace(k, v)
        write(self.root / page, out)
        self.pages.append(page)

    # --- sources ---------------------------------------------------------------------

    def markdown_page(self, md_path, page, placeholder_title, description="", meta_first=False, transform=None):
        md_path = Path(md_path)
        if md_path.exists():
            md = public.filter_markdown(read(md_path))   # the public filter: no When/Where rows, no forge links
            if transform:
                md = transform(md)
            body = pandoc(md)
            desc = description or first_sentence(body)
            self.render(page, placeholder_title, body, desc, meta_first=meta_first)
        else:
            self.notes.append("%s: source %s does not exist; placeholder written" % (page, md_path))
            body = "<h1>%s</h1>\n<p>This page has not been written yet.</p>\n" % html.escape(placeholder_title)
            self.render(page, placeholder_title, body, description)

    def reference_page(self):
        ref_dir = SITE / "reference"
        frags = sorted(p for p in ref_dir.glob("*.html"))
        page = "reference.html"
        if not frags:
            self.notes.append("reference.html: no fragments in site/reference/; placeholder written")
            body = "<h1>Reference</h1>\n<p>The function reference has not been generated yet (site/reference.R).</p>\n"
            self.render(page, "Reference", body)
            return
        index = None
        if (ref_dir / "index.json").exists():
            index = json.loads(read(ref_dir / "index.json"))
        if index is None:
            exports = set(re.findall(r"^export\(([^)]+)\)", read(self.root / "NAMESPACE"), flags=re.M))
            index = []
            for p in frags:
                topic = p.stem
                index.append({"topic": topic, "title": "", "kind": "export" if topic in exports else "internal", "file": str(p)})
            index.sort(key=lambda e: (e["kind"] != "export", e["topic"]))
        n_export = sum(1 for e in index if e.get("kind") == "export")
        parts = ["<h1>Reference</h1>"]
        parts.append(
            "<p>The documented functions of <code>emphasis</code> %s: the %d exported functions first, then the methods and the internal functions that the package documents.</p>"
            % (html.escape(self.version), n_export)
        )
        groups = [("export", "Exported functions"), ("S3method", "Methods"), ("internal", "Internal functions")]
        for kind, heading in groups:
            entries = [e for e in index if e.get("kind", "internal") == kind]
            if not entries:
                continue
            parts.append('<h2 id="%s">%s</h2>' % (kind, html.escape(heading)))
            for e in entries:
                topic = e["topic"]
                f = ref_dir / (topic + ".html")
                if not f.exists():
                    self.notes.append("reference.html: index.json names %s but site/reference/%s.html is missing" % (topic, topic))
                    continue
                frag = read(f)
                # Cross-references between topics are href="#<topic>" (site/reference.R); a link
                # into another package's help (stats.html#topic+optim) has no target here and
                # is reduced to its text.
                frag = re.sub(r'href="\.\./\.\./emphasis/help/([^"]+)"', r'href="#\1"', frag)
                frag = re.sub(r'<a href="(?!#)[^"]*">(.*?)</a>', r"\1", frag, flags=re.S)
                frag = re.sub(r"<h3>", '<h4>', frag)
                frag = re.sub(r"</h3>", "</h4>", frag)
                title = e.get("title") or ""
                all_aliases = e.get("aliases", [])
                # The Rd topic of an internal function is "dot-name"; the function is ".name".
                name = all_aliases[0] if all_aliases else topic
                aliases = [a for a in all_aliases if a != name]
                head = '<h3 id="%s"><code>%s</code>%s</h3>' % (
                    html.escape(topic), html.escape(name),
                    (" <span class=\"ref-title\">" + html.escape(title) + "</span>") if title else "")
                if aliases:
                    head += '<p class="ref-aliases">Also: %s</p>' % ", ".join("<code>%s</code>" % html.escape(a) for a in aliases)
                parts.append('<section class="ref">\n%s\n%s\n</section>' % (head, frag))
        self.render(page, "Reference", "\n".join(parts), "The functions of the emphasis package.",
                    numbered_allowed=False, run_in_toc=True)

    def experiments_index(self):
        home = self.wiki / "Home.md"
        if not home.exists():
            die("wiki has no Home.md")
        md = read(home)
        md = re.sub(r"^# .*$", "# Experiments", md, count=1, flags=re.M)
        # "The map" groups its tables under h3; they become the sections of this page.
        md = re.sub(r"^## The map\s*$", "", md, flags=re.M)
        md = re.sub(r"^### ", "## ", md, flags=re.M)
        body = pandoc(md)
        self.render("experiments/index.html", "Experiments", body,
                    "The experiments behind the emphasis package: the question, the status and the result of each, grouped by the question they answer.")

    def copy_assets(self):
        css_src = SITE / "style.css"
        if not css_src.exists():
            die("site/style.css is missing")
        write(self.root / "assets" / "site.css", read(css_src))
        if self.wiki:
            img_dir = self.wiki / "images"
            dest = self.root / "assets" / "wiki"
            dest.mkdir(parents=True, exist_ok=True)
            n = 0
            for p in sorted(img_dir.glob("*")):
                if p.is_file() and not p.name.startswith("."):
                    shutil.copy2(p, dest / p.name)
                    n += 1
            print("copied %d wiki images" % n)


def first_sentence(body):
    m = re.search(r"<p>(.*?)</p>", body, flags=re.S)
    if not m:
        return ""
    s = strip_tags(m.group(1))
    s = re.sub(r"\s+", " ", s)
    s = re.split(r"(?<=[.!?])\s", s)[0]
    return s[:300]


# ----------------------------------------------------------------------------- checker

ATTR_RE = re.compile(r'(?:href|src)="([^"]*)"')


def check_links(root, pages):
    """Every internal href/src in `pages` must resolve; fragments must exist on the target."""
    root = Path(root)
    ids = {}

    def page_ids(p):
        if p not in ids:
            try:
                t = read(root / p)
            except Exception:
                ids[p] = set()
                return ids[p]
            ids[p] = set(re.findall(r'\sid="([^"]+)"', t)) | set(re.findall(r'\sname="([^"]+)"', t))
        return ids[p]

    broken = []
    for page in pages:
        text = read(root / page)
        text = re.sub(r"<script.*?</script>", "", text, flags=re.S)
        for url in ATTR_RE.findall(text):
            if not url or re.match(r"^[a-z][a-z0-9+.-]*:", url) or url.startswith("//"):
                continue
            path, _, frag = url.partition("#")
            if path == "":
                target = page
            else:
                target = os.path.normpath(os.path.join(os.path.dirname(page), path))
            if not (root / target).exists():
                broken.append("%s -> %s (missing file)" % (page, url))
                continue
            if frag and target.endswith(".html") and frag not in page_ids(target):
                broken.append("%s -> %s (no such anchor)" % (page, url))
    return broken


# ----------------------------------------------------------------------------- main

def main():
    ap = argparse.ArgumentParser(description="build the emphasis site")
    ap.add_argument("--wiki", required=True, help="path to the experiments wiki checkout")
    ap.add_argument("--out", default=str(ROOT_DEFAULT), help="repository root to write into")
    args = ap.parse_args()

    root = Path(args.out).resolve()
    wiki = Path(args.wiki).resolve()
    if not wiki.is_dir():
        die("wiki path %s is not a directory" % wiki)
    for tool in ("pandoc",):
        if shutil.which(tool) is None:
            die(tool + " is not on PATH")
    template_path = SITE / "template.html"
    if not template_path.exists():
        die("site/template.html is missing")
    template = read(template_path)
    version = "0.0"
    m = re.search(r"^Version:\s*(\S+)", read(root / "DESCRIPTION"), flags=re.M)
    if m:
        version = m.group(1)

    site = Site(root, wiki, template, "version " + version)
    site.copy_assets()

    # Root pages from site/src.
    site.markdown_page(SITE / "src" / "index.md", "index.html", "Overview",
                       "The emphasis R package: diversification models with lineage-level covariates fitted to reconstructed phylogenies by Monte Carlo EM.")
    site.markdown_page(SITE / "src" / "usage.md", "usage.html", "Usage",
                       "Worked examples with the emphasis package and their output.")
    site.reference_page()

    # Wiki overview pages.
    for name, page in WIKI_PAGES.items():
        if name == "Home":
            continue
        site.markdown_page(wiki / (name + ".md"), page, name.replace("-", " "))

    # The experiments.
    site.experiments_index()
    exp = sorted((p for p in wiki.glob("*.md") if EXPERIMENT_RE.match(p.stem)),
                 key=lambda p: (p.stem != "S", int(p.stem[1:]) if p.stem != "S" else 0))
    for p in exp:
        site.markdown_page(p, "experiments/" + p.stem + ".html", p.stem, meta_first=True)
    print("wrote %d pages" % len(site.pages))

    for note in site.notes:
        print("note: " + note)

    broken = check_links(root, site.pages)
    if broken:
        print("broken links or images (%d):" % len(broken))
        for b in broken:
            print("  " + b)
        sys.exit(1)
    print("all internal links and images resolve on %d generated pages" % len(site.pages))
    # The public audit: no machine name, tailnet address, local path or library on any page.
    residue = []
    for p in site.pages:
        residue += public.audit(read(root / p), p)
    if residue:
        print("private tokens on public pages (%d):" % len(residue))
        for r in residue:
            print("  " + r)
        sys.exit(1)
    print("public audit: clean")
    # The kept pages are not edited by the build; a link from them that no longer resolves is
    # reported as a warning, since the fix belongs in the page they point to.
    kept = [k for k in KEPT if (root / k).exists()]
    for b in check_links(root, kept):
        print("warning: " + b)


if __name__ == "__main__":
    main()
