from __future__ import annotations

import html
import shutil
import zipfile
from pathlib import Path
from PIL import Image


ROOT = Path("/Users/sm3466/YSE Dropbox/Sparkle Malone/Research/FluxGradient/lterwg-flux-gradient-methane/workflows")
OUT = ROOT / "OUTPUT" / "NEON_CH4_Methods_Results_Presentation.pptx"

SLIDE_W = 12192000
SLIDE_H = 6858000
EMU_PER_IN = 914400


def emu(inches: float) -> int:
    return int(inches * EMU_PER_IN)


def esc(text: str) -> str:
    return html.escape(text, quote=True)


def rels_xml(rels: list[tuple[str, str, str]]) -> str:
    body = "\n".join(
        f'<Relationship Id="{rid}" Type="{typ}" Target="{esc(target)}"/>'
        for rid, typ, target in rels
    )
    return f'''<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">
{body}
</Relationships>'''


def text_body(lines: list[str], font_size: int, color: str = "1F2933", bold_first: bool = False) -> str:
    paragraphs = []
    for i, line in enumerate(lines):
        if line == "":
            paragraphs.append("<a:p/>")
            continue
        bold = ' b="1"' if (bold_first and i == 0) else ""
        paragraphs.append(
            f'''<a:p>
  <a:r><a:rPr lang="en-US" sz="{font_size * 100}"{bold}><a:solidFill><a:srgbClr val="{color}"/></a:solidFill></a:rPr><a:t>{esc(line)}</a:t></a:r>
  <a:endParaRPr lang="en-US" sz="{font_size * 100}"/>
</a:p>'''
        )
    return "\n".join(paragraphs)


def bullet_body(items: list[str], font_size: int, color: str = "1F2933") -> str:
    paragraphs = []
    for item in items:
        paragraphs.append(
            f'''<a:p>
  <a:pPr marL="228600" indent="-142875"><a:buChar char="•"/></a:pPr>
  <a:r><a:rPr lang="en-US" sz="{font_size * 100}"><a:solidFill><a:srgbClr val="{color}"/></a:solidFill></a:rPr><a:t>{esc(item)}</a:t></a:r>
  <a:endParaRPr lang="en-US" sz="{font_size * 100}"/>
</a:p>'''
        )
    return "\n".join(paragraphs)


def textbox_xml(shape_id: int, x: int, y: int, w: int, h: int, body: str, name: str = "Text") -> str:
    return f'''<p:sp>
<p:nvSpPr><p:cNvPr id="{shape_id}" name="{esc(name)}"/><p:cNvSpPr txBox="1"/><p:nvPr/></p:nvSpPr>
<p:spPr><a:xfrm><a:off x="{x}" y="{y}"/><a:ext cx="{w}" cy="{h}"/></a:xfrm><a:prstGeom prst="rect"><a:avLst/></a:prstGeom><a:noFill/><a:ln><a:noFill/></a:ln></p:spPr>
<p:txBody><a:bodyPr wrap="square" anchor="t"><a:spAutoFit/></a:bodyPr><a:lstStyle/>{body}</p:txBody>
</p:sp>'''


def rect_xml(shape_id: int, x: int, y: int, w: int, h: int, fill: str, line: str | None = None, name: str = "Block") -> str:
    ln = f'<a:ln w="9525"><a:solidFill><a:srgbClr val="{line}"/></a:solidFill></a:ln>' if line else '<a:ln><a:noFill/></a:ln>'
    return f'''<p:sp>
<p:nvSpPr><p:cNvPr id="{shape_id}" name="{esc(name)}"/><p:cNvSpPr/><p:nvPr/></p:nvSpPr>
<p:spPr><a:xfrm><a:off x="{x}" y="{y}"/><a:ext cx="{w}" cy="{h}"/></a:xfrm><a:prstGeom prst="rect"><a:avLst/></a:prstGeom><a:solidFill><a:srgbClr val="{fill}"/></a:solidFill>{ln}</p:spPr>
</p:sp>'''


def picture_xml(shape_id: int, rid: str, x: int, y: int, w: int, h: int, name: str = "Picture") -> str:
    return f'''<p:pic>
<p:nvPicPr><p:cNvPr id="{shape_id}" name="{esc(name)}"/><p:cNvPicPr><a:picLocks noChangeAspect="1"/></p:cNvPicPr><p:nvPr/></p:nvPicPr>
<p:blipFill><a:blip r:embed="{rid}"/><a:stretch><a:fillRect/></a:stretch></p:blipFill>
<p:spPr><a:xfrm><a:off x="{x}" y="{y}"/><a:ext cx="{w}" cy="{h}"/></a:xfrm><a:prstGeom prst="rect"><a:avLst/></a:prstGeom></p:spPr>
</p:pic>'''


def image_fit(path: Path, x_in: float, y_in: float, w_in: float, h_in: float) -> tuple[int, int, int, int]:
    with Image.open(path) as im:
        iw, ih = im.size
    box_w, box_h = emu(w_in), emu(h_in)
    scale = min(box_w / iw, box_h / ih)
    w, h = int(iw * scale), int(ih * scale)
    x = emu(x_in) + int((box_w - w) / 2)
    y = emu(y_in) + int((box_h - h) / 2)
    return x, y, w, h


class Slide:
    def __init__(self, title: str, subtitle: str | None = None):
        self.title = title
        self.subtitle = subtitle
        self.parts: list[str] = []
        self.images: list[tuple[Path, str]] = []
        self.next_id = 3

    def add_title(self):
        self.parts.append(rect_xml(2, 0, 0, SLIDE_W, emu(0.17), "0B3D3B"))
        self.parts.append(textbox_xml(self.next_id, emu(0.45), emu(0.27), emu(12.4), emu(0.55), text_body([self.title], 25, "0B3D3B", True), "Title"))
        self.next_id += 1
        if self.subtitle:
            self.parts.append(textbox_xml(self.next_id, emu(0.47), emu(0.82), emu(12.0), emu(0.35), text_body([self.subtitle], 10, "55606A"), "Subtitle"))
            self.next_id += 1

    def add_text(self, x, y, w, h, lines, size=14, color="1F2933", bold_first=False):
        self.parts.append(textbox_xml(self.next_id, emu(x), emu(y), emu(w), emu(h), text_body(lines, size, color, bold_first), "Text"))
        self.next_id += 1

    def add_bullets(self, x, y, w, h, items, size=14):
        self.parts.append(textbox_xml(self.next_id, emu(x), emu(y), emu(w), emu(h), bullet_body(items, size), "Bullets"))
        self.next_id += 1

    def add_card(self, x, y, w, h, title, body, accent="0B3D3B"):
        self.parts.append(rect_xml(self.next_id, emu(x), emu(y), emu(w), emu(h), "F7FAFA", "D6DEE0", "Card"))
        self.next_id += 1
        self.parts.append(rect_xml(self.next_id, emu(x), emu(y), emu(0.07), emu(h), accent, None, "Accent"))
        self.next_id += 1
        self.add_text(x + 0.18, y + 0.14, w - 0.35, 0.28, [title], 12, accent, True)
        self.add_text(x + 0.18, y + 0.48, w - 0.35, h - 0.55, body, 10, "1F2933")

    def add_image(self, path: Path, x, y, w, h):
        rid = f"rIdImg{len(self.images) + 1}"
        self.images.append((path, rid))
        ix, iy, iw, ih = image_fit(path, x, y, w, h)
        self.parts.append(picture_xml(self.next_id, rid, ix, iy, iw, ih, path.name))
        self.next_id += 1

    def add_footer(self, text: str):
        self.parts.append(textbox_xml(self.next_id, emu(0.45), emu(7.12), emu(12.3), emu(0.22), text_body([text], 7, "6B7280"), "Footer"))
        self.next_id += 1

    def xml(self) -> str:
        self.add_title()
        return f'''<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<p:sld xmlns:a="http://schemas.openxmlformats.org/drawingml/2006/main" xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships" xmlns:p="http://schemas.openxmlformats.org/presentationml/2006/main">
<p:cSld><p:bg><p:bgPr><a:solidFill><a:srgbClr val="FFFFFF"/></a:solidFill><a:effectLst/></p:bgPr></p:bg><p:spTree>
<p:nvGrpSpPr><p:cNvPr id="1" name=""/><p:cNvGrpSpPr/><p:nvPr/></p:nvGrpSpPr>
<p:grpSpPr><a:xfrm><a:off x="0" y="0"/><a:ext cx="0" cy="0"/><a:chOff x="0" y="0"/><a:chExt cx="0" cy="0"/></a:xfrm></p:grpSpPr>
{''.join(self.parts)}
</p:spTree></p:cSld><p:clrMapOvr><a:masterClrMapping/></p:clrMapOvr></p:sld>'''


def make_slides() -> list[Slide]:
    fig1 = ROOT / "FIGURES" / "NEON_DIEL_driver_flux_synthesis.png"
    fig2 = ROOT / "FIGURES" / "NEON_FLUXNET_CH4_flux_comparison.png"
    fig3 = ROOT / "FIGURES" / "NEON_strong_sink_driver_comparison.png"

    slides: list[Slide] = []

    s = Slide("NEON CH4 source-sink behavior", "Methods and results from half-hourly total-flux analysis")
    s.add_text(0.6, 1.45, 6.4, 1.1, ["Central claim", "NEON CH4 behavior is better interpreted as source-state switching plus separate emission and uptake magnitudes than as one continuous net-flux response."], 18, "0B3D3B", True)
    s.add_card(7.35, 1.25, 2.0, 1.2, "46", ["sites classified"], "0B3D3B")
    s.add_card(9.55, 1.25, 2.0, 1.2, "3", ["behavior classes"], "0B3D3B")
    s.add_card(7.35, 2.75, 2.0, 1.2, "9.1%", ["source-state deviance explained"], "0B3D3B")
    s.add_card(9.55, 2.75, 2.0, 1.2, "26-34%", ["magnitude-model deviance explained"], "0B3D3B")
    s.add_bullets(0.75, 3.45, 5.9, 1.8, [
        "Consistent sinks: 11 sites",
        "Fluctuating sites: 29 sites",
        "Consistent sources: 6 sites",
        "Absolute flux magnitudes remain provisional pending unit-conversion audit"
    ], 14)
    s.add_footer("Analysis outputs: OUTPUT/NEON_CH4_methods_results_draft.md")
    slides.append(s)

    s = Slide("Analysis workflow", "From half-hourly observations to site behavior classes and external benchmarks")
    steps = [
        ("Half-hour CH4 flux", ["flux_total = gradient flux + storage flux", "Converted to mg C m-2 30 min-1", "Positive = source; negative = sink"]),
        ("Site behavior", ["Monthly site summaries classify consistent sinks, fluctuating sites, sources"]),
        ("Two-part GAMs", ["Source probability", "Positive emission magnitude", "Sink uptake magnitude"]),
        ("Context", ["FLUXNET-CH4 ecosystem classes", "MeMo process-based upland uptake"]),
    ]
    xs = [0.55, 3.65, 6.75, 9.85]
    for x, (title, body) in zip(xs, steps):
        s.add_card(x, 1.55, 2.55, 3.45, title, body, "0B3D3B")
    s.add_text(0.8, 5.55, 11.6, 0.55, ["Key design choice: preserve half-hourly dynamics for diel/source-state models, but interpret site-level drivers using site summaries to avoid pseudoreplication."], 16, "1F2933")
    slides.append(s)

    s = Slide("Methods: source-state and magnitude models", "The analysis separates whether a site emits from how large emissions or uptake are")
    s.add_card(0.65, 1.25, 3.7, 2.1, "Source-state model", ["Binomial GAM", "Response: CH4_mgC_30min > 0", "Interprets probability of positive flux"], "0B3D3B")
    s.add_card(4.8, 1.25, 3.7, 2.1, "Emission magnitude", ["Gaussian GAM on log positive flux", "Conditional on CH4 flux > 0"], "1B6F6A")
    s.add_card(8.95, 1.25, 3.7, 2.1, "Uptake magnitude", ["Gaussian GAM on log absolute uptake", "Conditional on CH4 flux < 0"], "455A64")
    s.add_bullets(0.85, 4.0, 11.5, 1.65, [
        "Predictors: hour of day, sine/cosine diel terms, air temperature, VSWC, log(PAR), temperature x moisture, season, ecosystem type, CH4 behavior class, and random site effects.",
        "Apparent Q10 estimated within site-season-regime windows with at least 100 observations and at least 5 C temperature range.",
        "Implementation: R scripts NEON.DIEL.Analysis2.R and flow.30min.analysis.R."
    ], 13)
    slides.append(s)

    s = Slide("Result 1: source-state switching is the useful frame", "Continuous net flux alone obscures the structure in the data")
    s.add_card(0.65, 1.35, 3.2, 1.3, "Continuous flux model", ["Adjusted R2 = 0.003", "Deviance explained = 0.003"], "9CA3AF")
    s.add_card(4.05, 1.35, 3.2, 1.3, "Source-state model", ["Deviance explained = 9.09%"], "0B3D3B")
    s.add_card(7.45, 1.35, 2.35, 1.3, "Emissions", ["Deviance explained = 26%"], "1B6F6A")
    s.add_card(10.0, 1.35, 2.35, 1.3, "Uptake", ["Deviance explained = 33.6%"], "455A64")
    s.add_bullets(0.85, 3.35, 11.2, 1.7, [
        "Consistent sources are not just larger emitters; they enter positive-flux state more often.",
        "Consistent sinks stay in uptake-favorable states more often.",
        "Magnitude models become interpretable once positive emissions and sink uptake are separated."
    ], 15)
    s.add_footer("Figure support: NEON_DIEL_driver_flux_synthesis.png")
    slides.append(s)

    s = Slide("Figure 1: diel behavior and driver synthesis", "Source probability, magnitude, and Q10 tell different pieces of the CH4 story")
    s.add_image(fig1, 0.55, 1.15, 12.25, 5.65)
    s.add_footer("Figure file: FIGURES/NEON_DIEL_driver_flux_synthesis.png")
    slides.append(s)

    s = Slide("Result 2: apparent Q10 is modest, not a class separator", "Temperature sensitivity exists, but it does not explain source/sink identity")
    s.add_card(0.65, 1.25, 3.1, 1.65, "Positive emission Q10", ["Sinks: 1.28", "Fluctuating: 1.34", "Sources: 1.40", "Kruskal-Wallis p = 0.348"], "1B6F6A")
    s.add_card(4.05, 1.25, 3.1, 1.65, "Sink uptake Q10", ["Sinks: 1.46", "Fluctuating: 1.27", "Sources: 1.36", "Kruskal-Wallis p = 0.125"], "455A64")
    s.add_card(7.45, 1.25, 4.55, 1.65, "Interpretation", ["Q10 is an apparent temperature response, not causal proof.", "Moisture, diffusion, and site state remain central."], "0B3D3B")
    s.add_bullets(0.85, 3.65, 11.3, 1.4, [
        "Usable Q10 windows: 158 positive-emission site-seasons and 159 sink-uptake site-seasons.",
        "Q10 values were generally in the 1.3-1.5 range.",
        "The stronger manuscript point is state switching, not temperature sensitivity."
    ], 15)
    slides.append(s)

    s = Slide("Figure 2: NEON flux magnitudes in context", "NEON site medians are near zero relative to FLUXNET-CH4 and MeMo benchmarks")
    s.add_image(fig2, 0.45, 1.05, 12.45, 5.75)
    s.add_footer("References: Delwiche et al. 2021; Murguia-Flores et al. 2018")
    slides.append(s)

    s = Slide("Result 3: external comparisons expose a unit-conversion caveat", "Relative NEON classes are useful; absolute magnitudes need audit")
    s.add_bullets(0.75, 1.25, 5.8, 2.0, [
        "FLUXNET-CH4 ecosystem means are 7.95-127 mg C m-2 d-1.",
        "MeMo upland uptake estimates are about -1.24 to -0.216 mg C m-2 d-1.",
        "NEON source, fluctuating, and strong-sink medians are close to zero."
    ], 15)
    s.add_card(7.1, 1.25, 4.8, 1.1, "NEON source median", ["0.00262 mg C m-2 d-1"], "00008B")
    s.add_card(7.1, 2.55, 4.8, 1.1, "NEON fluctuating median", ["-0.000796 mg C m-2 d-1"], "666666")
    s.add_card(7.1, 3.85, 4.8, 1.1, "NEON strong sink median", ["-0.00635 mg C m-2 d-1"], "CC0000")
    s.add_text(0.75, 5.45, 11.4, 0.7, ["Caveat: these comparisons should be read as a diagnostic red flag for flux-unit verification before making claims about absolute NEON CH4 rates."], 16, "7A1F1F", True)
    slides.append(s)

    s = Slide("Figure 3: why are some sites stronger sinks?", "Strong sinks point toward moisture variability and soil physical controls")
    s.add_image(fig3, 0.55, 1.05, 12.25, 5.75)
    s.add_footer("Figure file: FIGURES/NEON_strong_sink_driver_comparison.png")
    slides.append(s)

    s = Slide("Interpretation and next steps", "A stronger story is emerging, with one necessary audit before publication")
    s.add_card(0.7, 1.25, 3.7, 2.1, "Main ecological result", ["Site behavior is mostly about source-state frequency: how often sites switch into positive CH4 flux."], "0B3D3B")
    s.add_card(4.8, 1.25, 3.7, 2.1, "Strong sink signal", ["Higher VSWC variance, lower clay, lower bulk density, higher sand: consistent with aeration/diffusion control."], "CC0000")
    s.add_card(8.9, 1.25, 3.7, 2.1, "Publication caveat", ["Absolute flux rates are tiny compared with literature; verify raw-unit conversion before rate claims."], "7A1F1F")
    s.add_bullets(0.9, 4.1, 11.2, 1.6, [
        "Keep source probability and magnitude as separate endpoints.",
        "Use Q10 as a secondary descriptive metric, not the central explanation.",
        "After unit audit, update literature comparison and manuscript text."
    ], 15)
    slides.append(s)

    s = Slide("References", "Primary references used in methods and comparison figures")
    s.add_bullets(0.7, 1.2, 12.0, 4.9, [
        "Delwiche et al. (2021). FLUXNET-CH4: a global, multi-ecosystem dataset and analysis of methane seasonality from freshwater wetlands. Earth System Science Data, 13, 3607-3689.",
        "Murguia-Flores et al. (2018). Soil Methanotrophy Model (MeMo v1.0): a process-based model to quantify global uptake of atmospheric methane by soil. Geoscientific Model Development, 11, 2009-2032.",
        "Wood (2017). Generalized Additive Models: An Introduction with R. Chapman and Hall/CRC.",
        "Repository outputs: NEON_CH4_methods_results_draft.md; NEON_DIEL_Analysis2_results.md; NEON_strong_sink_driver_comparison_results.md."
    ], 12)
    slides.append(s)

    return slides


def content_types(n_slides: int, image_names: list[str]) -> str:
    slide_overrides = "\n".join(
        f'<Override PartName="/ppt/slides/slide{i}.xml" ContentType="application/vnd.openxmlformats-officedocument.presentationml.slide+xml"/>'
        for i in range(1, n_slides + 1)
    )
    image_defaults = '<Default Extension="png" ContentType="image/png"/>'
    return f'''<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">
<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>
<Default Extension="xml" ContentType="application/xml"/>
{image_defaults}
<Override PartName="/ppt/presentation.xml" ContentType="application/vnd.openxmlformats-officedocument.presentationml.presentation.main+xml"/>
<Override PartName="/ppt/slideMasters/slideMaster1.xml" ContentType="application/vnd.openxmlformats-officedocument.presentationml.slideMaster+xml"/>
<Override PartName="/ppt/slideLayouts/slideLayout1.xml" ContentType="application/vnd.openxmlformats-officedocument.presentationml.slideLayout+xml"/>
<Override PartName="/ppt/theme/theme1.xml" ContentType="application/vnd.openxmlformats-officedocument.theme+xml"/>
{slide_overrides}
<Override PartName="/docProps/core.xml" ContentType="application/vnd.openxmlformats-package.core-properties+xml"/>
<Override PartName="/docProps/app.xml" ContentType="application/vnd.openxmlformats-officedocument.extended-properties+xml"/>
</Types>'''


def presentation_xml(n_slides: int) -> str:
    sld_ids = "\n".join(f'<p:sldId id="{255 + i}" r:id="rId{i}"/>' for i in range(1, n_slides + 1))
    return f'''<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<p:presentation xmlns:a="http://schemas.openxmlformats.org/drawingml/2006/main" xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships" xmlns:p="http://schemas.openxmlformats.org/presentationml/2006/main">
<p:sldMasterIdLst><p:sldMasterId id="2147483648" r:id="rId{n_slides + 1}"/></p:sldMasterIdLst>
<p:sldIdLst>{sld_ids}</p:sldIdLst>
<p:sldSz cx="{SLIDE_W}" cy="{SLIDE_H}" type="screen16x9"/>
<p:notesSz cx="6858000" cy="9144000"/>
<p:defaultTextStyle/>
</p:presentation>'''


def minimal_master() -> str:
    return '''<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<p:sldMaster xmlns:a="http://schemas.openxmlformats.org/drawingml/2006/main" xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships" xmlns:p="http://schemas.openxmlformats.org/presentationml/2006/main">
<p:cSld><p:bg><p:bgPr><a:solidFill><a:srgbClr val="FFFFFF"/></a:solidFill><a:effectLst/></p:bgPr></p:bg><p:spTree><p:nvGrpSpPr><p:cNvPr id="1" name=""/><p:cNvGrpSpPr/><p:nvPr/></p:nvGrpSpPr><p:grpSpPr><a:xfrm><a:off x="0" y="0"/><a:ext cx="0" cy="0"/><a:chOff x="0" y="0"/><a:chExt cx="0" cy="0"/></a:xfrm></p:grpSpPr></p:spTree></p:cSld>
<p:clrMap bg1="lt1" tx1="dk1" bg2="lt2" tx2="dk2" accent1="accent1" accent2="accent2" accent3="accent3" accent4="accent4" accent5="accent5" accent6="accent6" hlink="hlink" folHlink="folHlink"/>
<p:sldLayoutIdLst><p:sldLayoutId id="2147483649" r:id="rId1"/></p:sldLayoutIdLst>
<p:txStyles><p:titleStyle/><p:bodyStyle/><p:otherStyle/></p:txStyles>
</p:sldMaster>'''


def minimal_layout() -> str:
    return '''<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<p:sldLayout xmlns:a="http://schemas.openxmlformats.org/drawingml/2006/main" xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships" xmlns:p="http://schemas.openxmlformats.org/presentationml/2006/main" type="blank" preserve="1">
<p:cSld name="Blank"><p:spTree><p:nvGrpSpPr><p:cNvPr id="1" name=""/><p:cNvGrpSpPr/><p:nvPr/></p:nvGrpSpPr><p:grpSpPr><a:xfrm><a:off x="0" y="0"/><a:ext cx="0" cy="0"/><a:chOff x="0" y="0"/><a:chExt cx="0" cy="0"/></a:xfrm></p:grpSpPr></p:spTree></p:cSld>
<p:clrMapOvr><a:masterClrMapping/></p:clrMapOvr>
</p:sldLayout>'''


def theme() -> str:
    return '''<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<a:theme xmlns:a="http://schemas.openxmlformats.org/drawingml/2006/main" name="NEON CH4">
<a:themeElements><a:clrScheme name="NEON"><a:dk1><a:srgbClr val="0B3D3B"/></a:dk1><a:lt1><a:srgbClr val="FFFFFF"/></a:lt1><a:dk2><a:srgbClr val="1F2933"/></a:dk2><a:lt2><a:srgbClr val="F7FAFA"/></a:lt2><a:accent1><a:srgbClr val="0B3D3B"/></a:accent1><a:accent2><a:srgbClr val="CC0000"/></a:accent2><a:accent3><a:srgbClr val="00008B"/></a:accent3><a:accent4><a:srgbClr val="666666"/></a:accent4><a:accent5><a:srgbClr val="1B6F6A"/></a:accent5><a:accent6><a:srgbClr val="455A64"/></a:accent6><a:hlink><a:srgbClr val="1B6F6A"/></a:hlink><a:folHlink><a:srgbClr val="455A64"/></a:folHlink></a:clrScheme><a:fontScheme name="Aptos"><a:majorFont><a:latin typeface="Aptos Display"/></a:majorFont><a:minorFont><a:latin typeface="Aptos"/></a:minorFont></a:fontScheme><a:fmtScheme name="Clean"><a:fillStyleLst><a:solidFill><a:schemeClr val="phClr"/></a:solidFill></a:fillStyleLst><a:lnStyleLst><a:ln w="9525"><a:solidFill><a:schemeClr val="phClr"/></a:solidFill></a:ln></a:lnStyleLst><a:effectStyleLst><a:effectStyle><a:effectLst/></a:effectStyle></a:effectStyleLst><a:bgFillStyleLst><a:solidFill><a:schemeClr val="phClr"/></a:solidFill></a:bgFillStyleLst></a:fmtScheme></a:themeElements><a:objectDefaults/><a:extraClrSchemeLst/></a:theme>'''


def write_pptx():
    slides = make_slides()
    media_map: dict[Path, str] = {}
    for s in slides:
        for path, _ in s.images:
            if path not in media_map:
                media_map[path] = f"image{len(media_map) + 1}.png"

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(OUT, "w", compression=zipfile.ZIP_DEFLATED) as z:
        z.writestr("[Content_Types].xml", content_types(len(slides), list(media_map.values())))
        z.writestr("_rels/.rels", rels_xml([
            ("rId1", "http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument", "ppt/presentation.xml"),
            ("rId2", "http://schemas.openxmlformats.org/package/2006/relationships/metadata/core-properties", "docProps/core.xml"),
            ("rId3", "http://schemas.openxmlformats.org/officeDocument/2006/relationships/extended-properties", "docProps/app.xml"),
        ]))
        z.writestr("docProps/core.xml", '''<?xml version="1.0" encoding="UTF-8" standalone="yes"?><cp:coreProperties xmlns:cp="http://schemas.openxmlformats.org/package/2006/metadata/core-properties" xmlns:dc="http://purl.org/dc/elements/1.1/" xmlns:dcterms="http://purl.org/dc/terms/" xmlns:dcmitype="http://purl.org/dc/dcmitype/" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"><dc:title>NEON CH4 Methods and Results</dc:title><dc:creator>Codex</dc:creator><cp:lastModifiedBy>Codex</cp:lastModifiedBy></cp:coreProperties>''')
        z.writestr("docProps/app.xml", f'''<?xml version="1.0" encoding="UTF-8" standalone="yes"?><Properties xmlns="http://schemas.openxmlformats.org/officeDocument/2006/extended-properties" xmlns:vt="http://schemas.openxmlformats.org/officeDocument/2006/docPropsVTypes"><Application>Microsoft PowerPoint</Application><PresentationFormat>On-screen Show (16:9)</PresentationFormat><Slides>{len(slides)}</Slides></Properties>''')
        z.writestr("ppt/presentation.xml", presentation_xml(len(slides)))
        pres_rels = [(f"rId{i}", "http://schemas.openxmlformats.org/officeDocument/2006/relationships/slide", f"slides/slide{i}.xml") for i in range(1, len(slides) + 1)]
        pres_rels.append((f"rId{len(slides) + 1}", "http://schemas.openxmlformats.org/officeDocument/2006/relationships/slideMaster", "slideMasters/slideMaster1.xml"))
        z.writestr("ppt/_rels/presentation.xml.rels", rels_xml(pres_rels))
        z.writestr("ppt/slideMasters/slideMaster1.xml", minimal_master())
        z.writestr("ppt/slideMasters/_rels/slideMaster1.xml.rels", rels_xml([
            ("rId1", "http://schemas.openxmlformats.org/officeDocument/2006/relationships/slideLayout", "../slideLayouts/slideLayout1.xml"),
            ("rId2", "http://schemas.openxmlformats.org/officeDocument/2006/relationships/theme", "../theme/theme1.xml"),
        ]))
        z.writestr("ppt/slideLayouts/slideLayout1.xml", minimal_layout())
        z.writestr("ppt/slideLayouts/_rels/slideLayout1.xml.rels", rels_xml([
            ("rId1", "http://schemas.openxmlformats.org/officeDocument/2006/relationships/slideMaster", "../slideMasters/slideMaster1.xml")
        ]))
        z.writestr("ppt/theme/theme1.xml", theme())

        for idx, s in enumerate(slides, 1):
            z.writestr(f"ppt/slides/slide{idx}.xml", s.xml())
            rels = [("rIdLayout", "http://schemas.openxmlformats.org/officeDocument/2006/relationships/slideLayout", "../slideLayouts/slideLayout1.xml")]
            for path, rid in s.images:
                rels.append((rid, "http://schemas.openxmlformats.org/officeDocument/2006/relationships/image", f"../media/{media_map[path]}"))
            z.writestr(f"ppt/slides/_rels/slide{idx}.xml.rels", rels_xml(rels))

        for path, media_name in media_map.items():
            z.write(path, f"ppt/media/{media_name}")


if __name__ == "__main__":
    write_pptx()
    print(f"Wrote {OUT}")
