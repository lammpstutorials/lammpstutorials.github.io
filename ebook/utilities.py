def detect_block(n, file_content):
    within_block = True
    words_in_block = []
    while within_block:
        line = file_content[n]
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
    return words_in_block

def block_to_sentence(lines):
    #sentences = []
    #sentence = ''
    #for line in lines:
    #    words = line.split()
    #    for word in words:
    #        if word == '\n':
    #            sentences.append(sentence)
    #            sentence = ''
    #        else:
    #            sentence += word + ' '
    #sentences.append(sentence[:-1])
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