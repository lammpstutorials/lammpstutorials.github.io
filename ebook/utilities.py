def detect_block(n, file_content, keep_line_break=False):
    within_block = True
    words_in_block = []
    all_jump = []
    while within_block:
        try:
            line = file_content[n]
            jump = 0
            for letter in line:
                if letter == ' ':
                    jump += 1
                else:
                    break
            all_jump.append(jump)
            n += 1
            spitted_line = line.split()
            if line == '':
                words_in_block.append('\n')
            elif line[:2] == '  ':
                within_block = True
                for word in spitted_line:
                    words_in_block.append(word)
            else:
                within_block = False
            if keep_line_break:
                words_in_block.append('\n')
        except:
            within_block = False
    print()
    print(all_jump)
    print()
    return words_in_block

def block_to_sentence(lines):
    sentences = []
    sentence = ''
    for word in lines:
        if word == '\n':
            if (sentence != '\n') & (sentence != ''):
                sentences.append(sentence[:-1])
            sentence = ''
        else:
            sentence += word + ' '
    if (sentence != '\n') & (sentence != ''):
        sentences.append(sentence[:-1])
    return sentences

def count_occurence(main_string, sub_string):
    count=0
    start_index=0
    positions = []
    for _ in range(len(main_string)):
        j = main_string.find(sub_string,start_index)
        if(j!=-1):
            start_index = j+1
            count+=1
            positions.append(j)
    return count, positions

def extract_link(RST, label):
    n = 0
    position_link = []
    for m, line in enumerate(RST.file_content):
        if label in line:
            n += 1
            position_link.append(m)
    assert n == 2
    position_link = position_link[1]
    for sub in RST.file_content[position_link+2].split():
        if 'href' in sub:
            link = sub[6:-1]
    
    for sub in RST.file_content[position_link+2].split('>'):
        if '</a' in sub:
            text = sub[:-3]
    return link, text

def fix_link(RST, line):
    """Deal with links"""
    if '|' in line:
        _, pos_bar = count_occurence(line, '|')
        links = []
        texts = []
        labels = []
        for ini, end in zip(pos_bar[::2], pos_bar[1::2]):
            label = line[ini:end+1]
            labels.append(label[1:-1])
            link, text = extract_link(RST, label)
            links.append(link)
            texts.append(text)
        l = 0
        sentence = line.split('|')
        new_line = ''
        for part in sentence:
            if part not in labels:
                new_line += part
                try:
                    replacement = '\href{'+links[l]+'}{'+texts[l]+'}'
                    new_line += replacement
                except:
                    pass                        
                l += 1
        return new_line
    else:
        return line

def fix_caption(lines):
    new_lines = []
    caption = None
    for line in lines:
        if 'caption' in line:
            caption = line.split('*')[1]+':'
        else:
            new_lines.append(line)
    return caption, new_lines


def fix_italic(line, caracter = '*', replace_with = [r'\textit{', '}'], replace_underscore=False):
    """Deal with special characters"""
    if caracter in line:
        _, pos_bar = count_occurence(line, caracter)
        italics = []
        for ini, end in zip(pos_bar[::2], pos_bar[1::2]):
            italics.append(line[ini:end+1][1:-1])
        l = 0
        sentence = line.split(caracter)
        new_line = ''
        for part in sentence:
            if part not in italics:
                new_line += part
                try:
                    italic = italics[l]
                    if replace_underscore:
                        italic = replace_special_character(italic, '_', '$\_$')
                    replacement = replace_with[0]+italic+replace_with[1]
                    new_line += replacement
                except:
                    pass                        
                l += 1
        return new_line
    else:
        return line
    
def fix_math(line):
    if ':math:' in line:
        _, pos_math = count_occurence(line, ':math:')
        _, pos_eq = count_occurence(line, "`")
        equations = []
        for ini, end in zip(pos_eq[::2], pos_eq[1::2]):
            equations.append(line[ini:end+1][1:-1])
        l = 0
        sentence = line.split("`")
        new_line = ''
        for part in sentence:
            if part not in equations:
                new_line += part
            else:
                try:
                    new_line += '$'+equations[l]+'$'             
                    l += 1
                except:
                    pass
        rests = new_line.split(':math:')
        new_line = ''
        for rest in rests:
            new_line += rest
        return new_line
    else:
        return line

def replace_special_character(line, to_replace, replace_with):
    if to_replace in line:
        sentence = line.split(to_replace)
        if len(line.split(to_replace)) == 2:
            new_line = sentence[0] + replace_with + sentence[1]
        else:
            print(len(line.split(to_replace)))          
            stop
        return new_line
    else:
        return line
