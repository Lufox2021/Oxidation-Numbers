# Rebuilds index.html from src/: inlines chem.js + org.js into template.html
from pathlib import Path
s = Path('src')
t = (s/'template.html').read_text()
c = (s/'chem.js').read_text() + (s/'org.js').read_text()
Path('index.html').write_text(t.replace('/*CHEM*/', c))
Path('tests/all.js').write_text(c)
print('index.html rebuilt')
