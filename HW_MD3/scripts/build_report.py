"""Render this project's Markdown report as an offline, printable HTML file.

Small standard-library renderer for the specific report syntax; no network assets.
"""
from pathlib import Path
import base64
import html
import re
import io
import matplotlib
matplotlib.use('Agg')
from matplotlib.mathtext import math_to_image

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / 'report/Cu表面能与低温平衡形貌.md'

def math_image(expression, display=False):
    buffer = io.BytesIO()
    math_to_image('$'+expression+'$', buffer, format='svg', dpi=160)
    encoded = base64.b64encode(buffer.getvalue()).decode()
    cls = 'display-math' if display else 'inline-math'
    return f'<img class="{cls}" alt="{html.escape(expression, quote=True)}" src="data:image/svg+xml;base64,{encoded}">'

def inline(text):
    formulas = []
    def capture(match):
        formulas.append(math_image(match[1]))
        return f'OVITOMATHTOKEN{len(formulas)-1}END'
    text = re.sub(r'(?<!\\)\$([^$\n]+)\$', capture, text)
    text = html.escape(text)
    text = re.sub(r'`([^`]+)`', r'<code>\1</code>', text)
    text = re.sub(r'\*\*([^*]+)\*\*', r'<strong>\1</strong>', text)
    text = re.sub(r'\*([^*]+)\*', r'<em>\1</em>', text)
    text = re.sub(r'\[([^\]]+)\]\(([^)]+)\)', r'<a href="\2">\1</a>', text)
    for i, formula in enumerate(formulas):
        text = text.replace(f'OVITOMATHTOKEN{i}END', formula)
    return text

lines = SOURCE.read_text().splitlines()
blocks, paragraph = [], []

def flush():
    if paragraph:
        blocks.append('<p>'+'<br>'.join(inline(p.strip()) for p in paragraph)+'</p>')
        paragraph.clear()

i = 0
while i < len(lines):
    line = lines[i]
    if not line.strip():
        flush()
    elif line.strip() == '$$':
        flush()
        formula = []
        i += 1
        while i < len(lines) and lines[i].strip() != '$$':
            formula.append(lines[i])
            i += 1
        blocks.append('<div class="equation">'+math_image(' '.join(formula), display=True)+'</div>')
    elif line.startswith('> '):
        flush()
        quote = []
        while i < len(lines) and lines[i].startswith('>'):
            quote.append(inline(lines[i].removeprefix('>').strip()))
            i += 1
        blocks.append('<blockquote>'+'<br>'.join(quote)+'</blockquote>')
        i -= 1
    elif line.startswith('#'):
        flush()
        level = len(line)-len(line.lstrip('#'))
        blocks.append(f'<h{level}>{inline(line[level:].strip())}</h{level}>')
    elif line.startswith('```'):
        flush()
        code = []
        i += 1
        while i < len(lines) and not lines[i].startswith('```'):
            code.append(lines[i])
            i += 1
        blocks.append('<pre><code>'+html.escape('\n'.join(code))+'</code></pre>')
    elif line.startswith('|'):
        flush()
        table = []
        while i < len(lines) and lines[i].startswith('|'):
            table.append([x.strip() for x in lines[i].strip('|').split('|')])
            i += 1
        if len(table)<2 or not all(re.fullmatch(r':?-+:?', x) for x in table[1]):
            raise ValueError('Malformed table')
        blocks.append('<table><thead><tr>'+''.join('<th>'+inline(x)+'</th>' for x in table[0])+'</tr></thead><tbody>')
        for row in table[2:]:
            if len(row)!=len(table[0]):
                raise ValueError('Table column mismatch')
            blocks.append('<tr>'+''.join('<td>'+inline(x)+'</td>' for x in row)+'</tr>')
        blocks.append('</tbody></table>')
        i -= 1
    elif line.startswith('!['):
        flush()
        match = re.fullmatch(r'!\[([^]]*)\]\(([^)]+)\)', line)
        if not match:
            raise ValueError('Malformed image')
        path = (SOURCE.parent / match[2]).resolve()
        data = base64.b64encode(path.read_bytes()).decode()
        blocks.append(f'<figure><img alt="{html.escape(match[1])}" src="data:image/png;base64,{data}"></figure>')
    elif re.match(r'^\d+\. ',line):
        flush()
        blocks.append('<p class="numbered">'+inline(line)+'</p>')
    else:
        paragraph.append(line)
    i += 1
flush()
css = '''
:root{color-scheme:light}*{box-sizing:border-box}body{margin:0;background:#eff2f5;color:#243344;
font-family:"Noto Sans CJK SC","Source Han Sans CN",sans-serif;font-size:15px;line-height:1.85}
main{max-width:980px;margin:32px auto;background:white;padding:56px 64px;box-shadow:0 8px 36px #24334412}
h1{font-size:29px;line-height:1.45;color:#173d61;margin:0 0 24px;border-bottom:3px solid #5083ab;padding-bottom:22px}
h2{font-size:21px;color:#173d61;margin:34px 0 14px}h3{font-size:17px;margin:24px 0 10px}
main>p:first-child{font-size:27px;line-height:1.5;border-bottom:3px solid #5083ab;padding-bottom:20px}
blockquote{margin:18px 0;padding:12px 18px;border-left:3px solid #5083ab;background:#f3f6fa}
.equation{text-align:center;margin:20px 0;break-inside:avoid}
img.display-math{display:inline-block;width:auto;max-width:100%;height:1.9em;object-fit:contain}
img.inline-math{display:inline-block;width:auto;height:1.1em;vertical-align:-.13em}
p{margin:12px 0;text-align:justify;overflow-wrap:anywhere}strong{color:#173d61}
table{width:100%;border-collapse:collapse;font-size:12px;line-height:1.6;margin:20px 0}
th{background:#eaf0f6;color:#173d61}th,td{padding:9px 8px;border-bottom:1px solid #d6e0e8;text-align:left}
tbody tr:nth-child(even){background:#f7f9fb}figure{margin:22px 0 8px}img{display:block;width:100%;height:auto}
code{font-family:monospace;font-size:.86em;background:#edf1f4;padding:1px 4px;overflow-wrap:anywhere}
pre{padding:14px;background:#edf1f4;white-space:pre-wrap;line-height:1.6}pre code{padding:0}
a{color:#28658c;text-decoration:none}a:hover{text-decoration:underline}
@media(max-width:760px){main{padding:28px 20px;margin:0}table{font-size:10px}h1{font-size:24px}}
@page{size:A4;margin:18mm 16mm}
@media print{body{background:white;font-size:10pt;line-height:1.65}main{margin:0;padding:0;max-width:none;box-shadow:none}
h1{font-size:20pt}h2{font-size:14pt}h3{font-size:11.5pt}h1,h2,h3{break-after:avoid}
figure,table,pre{break-inside:avoid}table{font-size:8pt}a{color:inherit}img{max-height:92mm;object-fit:contain}
p{orphans:3;widows:3}}
'''
output = '<!doctype html><html lang="zh-CN"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>FCC Cu 表面能与低温平衡形貌</title><style>'+css+'</style></head><body><main>'+''.join(blocks)+'</main></body></html>'
target = SOURCE.with_suffix('.html')
target.write_text(output)
print(target)
