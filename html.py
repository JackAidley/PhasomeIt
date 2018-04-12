'HTML helper class'

def makeLink(link, linkName=None) :
    if not linkName :
        name = link
    else :
        name = linkName
    return '<a href="'+link+'">'+name+'</a>'

def makeTag(tag, data, attributes=None) :
    attribText = ""
    if attributes :
        for attrib, value in attributes.items() :
            attribText += ' '+attrib+'="'+str(value)+'"'
    return "<"+tag+attribText+">"+str(data)+"</"+tag+">"

def sciFormat(number, format='.2e') :
    str = ('{0:'+format+'}').format(number)
    split = str.split('e')
    if len(split) == 1 :
        return split

    return split[0]+'&times;10<sup>'+split[1]+'</sup>'

class Styler :
    'Class for managing document style'

    def __init__(self) :
        self.styles = None

    def addStyle(self, styles) :
        if not self.styles :
            self.styles =[]
        self.styles = self.styles + styles

    def writeHeader(self, fp) :
        if self.styles :
            fp.write('<style media="screen" type="text/css">')
            for style in self.styles :
                fp.write(style)
            fp.write('</style>\n')

headerFunctions = {
# Function to add multiple functions to the windows load
'addLoadEvent' :
'''function addLoadEvent(func) {
  var oldonload = window.onload;
  if (typeof window.onload != 'function') {
    window.onload = func;
  } else {
    window.onload = function() {
      if (oldonload) {
        oldonload();
      }
      func();
    }
  }
}''',

'colouredText' :
'''function colouredText(canvas, col0, col1, text, x, y) {
  if (col0 != col1) {
    var oddText = "";
    var evenText = "";
    for( var i = 0; i < text.length; ++i) {
      if (i % 2 == 0) {
        oddText += text[i];
        evenText += " ";
      } else {
        oddText += " ";
        evenText += text[i];
      }
    }
    canvas.fillStyle = col0;
    canvas.fillText(oddText, x, y);
    canvas.fillStyle = col1;
    canvas.fillText(evenText, x, y);
  } else {
    canvas.fillStyle = col0;
    canvas.fillText(text, x, y);
  }
}''' }

class Doc :
    'Class for creating an HTML document'

    def __init__(self, path, title=None, styler=None, includeInHeader = None, sortTable=None) :
        self.fp = open(path, "w")
        self.styler=styler
        self.debugMode = False
        self.counter = 0
        self.sortTable=sortTable
        self.writeHeader(title, includeInHeader)

    def __enter__(self) :
        'Nothing'
        return self

    def __exit__(self, exception_type, exception_val, trace):
        self.fp.write("</body></html>")
        self.fp.close()

    def writeHeader(self, title, includeInHeader) :
        self.fp.write(
        """<!doctype html>
        <html lang="en">
        <head>
          <meta charset="utf-8">""")

        if title is not None:
            self.fp.write("<title>"+title+"</title>\n")

        if self.sortTable :
            self.fp.write('<script src="'+self.sortTable+'"></script>')

        if self.styler :
            self.styler.writeHeader(self.fp)

        if self.sortTable :
            self.fp.write('<style media="screen" type="text/css">')
            self.fp.write('table.sortable th:not(.sorttable_sorted):not(.sorttable_sorted_reverse):not(.sorttable_nosort):after { content: " \\25B4\\25BE"}')
            self.fp.write('</style>\n')

        if includeInHeader :
            with self.tag('script', attributes={'type' : 'text/javascript'}) :
                for fn in includeInHeader :
                    self.fp.write(headerFunctions[fn]+'\n')

        self.fp.write("</head><body>")

    def add(self,tag,data,attributes=None) :
        self.addRaw(makeTag(tag, data, attributes))

    def addExpandable(self, link, text) :
        with self.tag('script', attributes={'type' : 'text/javascript'}) :
            self.addRaw('function resizeIFrame'+str(self.counter)+'() {');
            self.addRaw('  var frame = window.parent.document.getElementById("insertFrame'+str(self.counter)+'");')
            self.addRaw('  frame.style.height = frame.contentWindow.document.body.scrollHeight + "px";')
            self.addRaw('}')
            self.addRaw('function expandSection'+str(self.counter)+'() {');
            self.addRaw('  document.getElementById("section'+str(self.counter)+'").innerHTML= "<iframe id=\'insertFrame'+str(self.counter)+'\' width=\'100%\' onload=\'resizeIFrame'+str(self.counter)+'()\' src=\''+link+'\' frameBorder=\'0\' scrolling=\'no\'><\/iframe>";')
            self.addRaw('};')
        with self.tag('div', attributes={'id' : 'section'+str(self.counter)}) :
            self.add('a', text, attributes={'href' : 'javascript:expandSection'+str(self.counter)+'()'})
        self.counter += 1

    # Naming here is bad, sorry
    # This one adds a hidden block of text
    # The above one imports a file into the hidden block
    def addCollapsed(self, hiddenText, text) :
        with self.tag('script', attributes={'type' : 'text/javascript'}) :
            self.addRaw('function resizeIFrame'+str(self.counter)+'() {');
            self.addRaw('  var frame = window.parent.document.getElementById("insertFrame'+str(self.counter)+'");')
            self.addRaw('  frame.style.height = frame.contentWindow.document.body.scrollHeight + "px";')
            self.addRaw('}')
            self.addRaw('function expandSection'+str(self.counter)+'() {');
            self.addRaw('  document.getElementById("section'+str(self.counter)+'").innerHTML= "'+hiddenText+'";')
            self.addRaw('};')
        with self.tag('div', attributes={'id' : 'section'+str(self.counter)}) :
            self.add('a', text, attributes={'href' : 'javascript:expandSection'+str(self.counter)+'()'})
        self.counter += 1


    def addRaw(self, text) :
        self.fp.write(text)
        if self.debugMode :
            self.fp.write('\n')

    def addRule(self) :
        self.addRaw('<hr />')

    def table(self, **kwargs) :
        return self.Table(self, **kwargs)

    def tag(self, tag, attributes=None) :
        return self.Tag(self, tag, attributes)

    # Table class
    class Table :
        def __init__(self, parent, sortable=None) :
            self.parent = parent
            self.firstRow = True
            self.bodyComplete = False
            if sortable is None :
                self.sortable = self.parent.sortTable
            else :
                self.sortable = sortable

        def __enter__(self) :
            if self.sortable :
                self.parent.fp.write('<table class="sortable">')
            else :
                self.parent.fp.write("<table>")
            return self

        def __exit__(self, exception_type, exception_val, trace):
            if not self.bodyComplete and not self.firstRow:
                self.parent.fp.write('</tbody>')
            self.parent.fp.write("</table>")

        def addHeaderRow(self, row) :
            self.parent.fp.write("<thead><tr>")
            for entry in row :
                self.parent.add("th", entry)
            self.parent.fp.write("</tr></thead>")

        def addRow(self, row, attributes = {}) :
            if self.firstRow :
                self.parent.fp.write("<tbody>")
                self.firstRow = False
            attribText=''
            for attrib, value in attributes.items() :
                attribText += ' '+attrib+'="'+str(value)+'"'
            self.parent.fp.write('<tr'+attribText+'>')
            for entry in row :
                self.parent.add('td', entry)
            self.parent.fp.write('</tr>')

        def addFooterRow(self, row) :
            self.bodyComplete=True
            if not self.firstRow :
                self.parent.fp.write('</tbody>')
            self.parent.fp.write("<tfoot><tr>")
            for entry in row :
                self.parent.add("td", entry)
            self.parent.fp.write("</tr></tfoot>")

        class Row :
            def __init__(self, parent, attributes={}) :
                self.parent = parent
                self.attributes = attributes

            def __enter__(self) :
                if self.parent.firstRow :
                    self.parent.parent.fp.write("<tbody>")
                    self.parent.firstRow = False
                attribText=''
                for attrib, value in self.attributes.items() :
                    attribText += ' '+attrib+'="'+str(value)+'"'
                self.parent.parent.fp.write('<tr'+attribText+'>')
                return self

            def __exit__(self, exception_type, exception_val, trace):
                self.parent.parent.fp.write('</tr>')

            def addCell(self, text, attributes = {}) :
                attribText=''
                for attrib, value in attributes.items() :
                    attribText += ' '+attrib+'="'+str(value)+'"'
                self.parent.parent.fp.write('<td'+attribText+'>'+text+'</td>')

        def row(self, **kwargs) :
            return self.Row(self, **kwargs)

    # Generic tag class
    class Tag :
        def __init__(self,  parent, tag, attributes) :
            self.tag = tag
            self.attributes = attributes
            self.parent = parent

        def __enter__(self) :
            attribText = ""
            if self.attributes :
                for attrib, value in self.attributes.items() :
                    attribText += ' '+attrib+'="'+str(value)+'"'
            self.parent.addRaw("<"+self.tag+attribText+">")
            return self

        def __exit__(self, exception_type, exception_val, trace) :
            self.parent.addRaw("</"+self.tag+">")