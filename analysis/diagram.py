"""Writes media/network.svg: the two merged pi-sections feeding the RL load (light card, dark text)."""
from common import MEDIA

W, H = 900, 300
els = []


def line(x1, y1, x2, y2):
    els.append(f'<line x1="{x1}" y1="{y1}" x2="{x2}" y2="{y2}" stroke="#222" stroke-width="2"/>')


def box(x, y, w, h, label, sub=""):
    els.append(f'<rect x="{x}" y="{y}" width="{w}" height="{h}" fill="#fff" stroke="#222" stroke-width="2"/>')
    els.append(f'<text x="{x + w / 2}" y="{y + h / 2 + 5}" text-anchor="middle" font-size="15">{label}</text>')
    if sub:
        els.append(f'<text x="{x + w / 2}" y="{y - 8}" text-anchor="middle" font-size="12" fill="#555">{sub}</text>')


def cap(x, label, sub):
    line(x, 150, x, 172); line(x - 14, 172, x + 14, 172); line(x - 14, 180, x + 14, 180); line(x, 180, x, 230)
    els.append(f'<text x="{x + 20}" y="180" font-size="14">{label}</text>')
    els.append(f'<text x="{x + 20}" y="196" font-size="11" fill="#555">{sub}</text>')


top = 150
line(60, top, 840, top); line(60, 230, 840, 230)
els.append('<circle cx="60" cy="190" r="20" fill="#fff" stroke="#222" stroke-width="2"/>')
els.append('<text x="60" y="195" text-anchor="middle" font-size="14">u1</text>')
els.append('<circle cx="840" cy="190" r="20" fill="#fff" stroke="#222" stroke-width="2"/>')
els.append('<text x="840" y="195" text-anchor="middle" font-size="14">u2</text>')
line(60, 150, 60, 170); line(60, 210, 60, 230); line(840, 150, 840, 170); line(840, 210, 840, 230)
# series elements drawn over the top rail
els.append(f'<rect x="150" y="{top-14}" width="60" height="28" fill="#fff" stroke="#222" stroke-width="2"/>')
els.append(f'<text x="180" y="{top-22}" text-anchor="middle" font-size="13">R1 = 30 Ω</text>')
els.append(f'<rect x="240" y="{top-14}" width="60" height="28" fill="#fff" stroke="#222" stroke-width="2"/>')
els.append(f'<text x="270" y="{top-22}" text-anchor="middle" font-size="13">L1 = 0.1 H</text>')
els.append(f'<rect x="600" y="{top-14}" width="60" height="28" fill="#fff" stroke="#222" stroke-width="2"/>')
els.append(f'<text x="630" y="{top-22}" text-anchor="middle" font-size="13">R2 = 30 Ω</text>')
els.append(f'<rect x="690" y="{top-14}" width="60" height="28" fill="#fff" stroke="#222" stroke-width="2"/>')
els.append(f'<text x="720" y="{top-22}" text-anchor="middle" font-size="13">L2 = 0.2 H</text>')
cap(110, "C1", "0.5 µF"); cap(340, "C1", "0.5 µF"); cap(560, "C2", "1.5 µF"); cap(790, "C2", "1.5 µF")
# load
els.append('<rect x="425" y="140" width="50" height="22" fill="#fff" stroke="#222" stroke-width="2"/>')
els.append('<text x="450" y="156" text-anchor="middle" font-size="12">RL</text>')
els.append('<rect x="425" y="192" width="50" height="22" fill="#fff" stroke="#222" stroke-width="2"/>')
els.append('<text x="450" y="208" text-anchor="middle" font-size="12">LL</text>')
line(450, 150, 450, 140); line(450, 162, 450, 192); line(450, 214, 450, 230)
els.append('<text x="450" y="118" text-anchor="middle" font-size="14" fill="#B03000">u_load (output, state x4)</text>')
els.append('<text x="450" y="268" text-anchor="middle" font-size="12" fill="#555">line 1 (100 km) | load RL = 100 Ω, LL = 1 H | line 2 (200 km)</text>')
els.append('<text x="450" y="30" text-anchor="middle" font-size="16" font-weight="bold">Two pi-model lines feeding a series RL load</text>')
els.append('<text x="450" y="52" text-anchor="middle" font-size="12" fill="#555">states: i1, i2, iL, u_load; inputs: u1, u2</text>')

svg = (f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" width="{W}" height="{H}" '
       f'font-family="DejaVu Sans, Arial, sans-serif" role="img" aria-label="Circuit of two pi-model transmission lines and an RL load">'
       f'<rect width="{W}" height="{H}" rx="12" fill="#f6f7f9"/>' + "".join(els) + "</svg>")
(MEDIA / "network.svg").write_text(svg)
print("wrote", MEDIA / "network.svg")
