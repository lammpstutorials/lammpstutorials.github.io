from utilities import is_in_verbatim, read_label, read_link

class FixDocument:
    """Fix Tex file."""
    def __init__(self, file_name, *args, **kwargs,):
        """Initialize"""
        super().__init__(*args, **kwargs)
        self.file_name = file_name

    def fix_document(self):
        """Improve the final document"""
        self.remove_space()
        keywords = [r'\end{lcverbatim}', r'\Large', r'\section{',
                    r'\subsection{', r'\end{wrapfigure}', r'\end{tcolorbox}']
        self.add_non_indent(keywords)
        self.convert_itemize()
        self.fix_label()
        self.fix_external_link()

    def add_non_indent(self, keywords):
        """Add nonindent command where appropriate"""
        initial_tex_file = self.import_tex_file()
        new_tex_file_name = []
        add_noindent = False
        for line in initial_tex_file:
            if line != '\n':
                prev_add_noindent = add_noindent
                add_noindent = False
                for word in keywords:
                    if word in line:
                        add_noindent = True
            if prev_add_noindent: 
                new_tex_file_name.append(r'\noindent '+line)
                add_noindent = False
                prev_add_noindent = False
            else:
                new_tex_file_name.append(line)
        self.write_file(new_tex_file_name)

    def remove_space(self):
        """Remove unecessary space between lines"""
        initial_tex_file = self.import_tex_file()
        # remove extra space
        new_tex_file_name = []
        consecutive_empty = 0
        for line in initial_tex_file:
            # remove double space
            if line == '\n':
                consecutive_empty += 1
            else:
                consecutive_empty = 0
            if consecutive_empty <= 1:
                new_tex_file_name.append(line)
        self.write_file(new_tex_file_name)

    def convert_itemize(self):
        """Convert list into latex item lists"""
        initial_tex_file = self.import_tex_file()
        in_verbatim = False
        add_start = False
        add_end = False
        within_item = False
        new_tex_file_name = []
        consecutive_item = 0
        for line in initial_tex_file:
            in_verbatim = is_in_verbatim(in_verbatim, line)
            add_start = False
            add_end = False
            if in_verbatim is False:
                if (line[0] == '-'):
                    line = '\item' + line[1:]
                    if consecutive_item == 0:
                        add_start = True
                    consecutive_item += 1
                else:
                    if (consecutive_item > 0) & (line[0] != ' '):
                        add_end = True
                    consecutive_item = 0
                if (add_start) & (within_item is False):
                    new_tex_file_name.append(r'\begin{itemize}'+'\n')
                    add_start = False
                    within_item = True
                elif (add_end) & (within_item):
                    new_tex_file_name.append(r'\end{itemize}'+'\n')
                    add_end = False
                    within_item = False
                    consecutive_item = 0
            new_tex_file_name.append(line)
        self.write_file(new_tex_file_name)

    def import_tex_file(self):
        """Import tex file"""
        f = open(self.file_name, "r") 
        initial_tex_file = []
        for line in f:
            initial_tex_file.append(line)
        f.close()
        return initial_tex_file

    def write_file(self, new_tex_file_name):
        """Write tex file"""
        f = open(self.file_name, "w") 
        for line in new_tex_file_name:
            f.write(line)
        f.close()

    def fix_label(self):
        """Fix label"""
        initial_tex_file = self.import_tex_file()
        # remove extra space
        new_tex_file_name = []
        for line in initial_tex_file:
            if ':ref:' in line:
                split = line.split(":ref:")
                if len(split) == 2:
                    label, rest = read_label(split[1])
                    if label is None:
                        print("WARNING, Wrong label", line)
                    new_line = split[0] + r'Tutorial\,\ref{' + label + '}' + rest[1]
                    new_tex_file_name.append(new_line)
                else:
                    print("WARNING, several label per line", line)
            else:
                new_tex_file_name.append(line)
        self.write_file(new_tex_file_name)

    def fix_external_link(self):
        """Fix external link"""
        initial_tex_file = self.import_tex_file()
        # remove extra space
        new_tex_file_name = []
        for line in initial_tex_file:
            if '../../inputs/' in line:
                rest, link = read_link(line)
                assert "\href" in rest[0], """Unexpected link"""
                actual_link = 'https://lammpstutorials.github.io/'+link[0].split('../')[-1]
                new_line = rest[0] + "{" +  actual_link + "}{"+ link[1] +  "}" + rest[2]
                #block = new_line.split("\href")
                #assert len(block) == 2
                #new_line = block[0] + r'Tutorial\,\href' + block[1]
                new_tex_file_name.append(new_line)
            else:
                new_tex_file_name.append(line)
        self.write_file(new_tex_file_name)      
