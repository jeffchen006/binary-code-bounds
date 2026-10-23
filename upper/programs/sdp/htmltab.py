from html.parser import HTMLParser
import re
class _P(HTMLParser):
    def __init__(s):
        super().__init__(); s.tables=[]; s.cur=None; s.row=None; s.cell=None; s.text=[]; s.insup=False
    def handle_starttag(s,tag,attrs):
        tag=tag.lower(); a={k.lower():v for k,v in attrs}
        if tag=='table': s.cur=[]; s.tables.append((''.join(s.text[-60:]),s.cur))
        elif tag=='tr':
            s.close_cell(); s.row=[]
            if s.cur is not None: s.cur.append(s.row)
        elif tag in('td','th'):
            s.close_cell(); s.cell=[tag,int(a.get('colspan',1) or 1),'','']
        elif tag=='sup': s.insup=True
    def close_cell(s):
        if s.cell is not None and s.row is not None: s.row.append(s.cell)
        s.cell=None
    def handle_endtag(s,tag):
        tag=tag.lower()
        if tag=='sup': s.insup=False
        if tag=='table': s.close_cell(); s.cur=None; s.row=None
    def handle_data(s,data):
        if s.cell is not None:
            if s.insup: s.cell[3]+=data
            else: s.cell[2]+=data
        else: s.text.append(data)
def tables(path):
    p=_P(); p.feed(open(path,encoding='utf-8',errors='replace').read()); p.close_cell()
    return [(re.sub(r'\s+',' ',c),t) for c,t in p.tables]
