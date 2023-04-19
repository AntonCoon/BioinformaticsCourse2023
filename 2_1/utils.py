class Read:
    def __init__(self, file):
        self.seq_id = int(read_nonempty_line(file)[1:])
        self.seq = read_nonempty_line(file)
        comment_list = read_nonempty_line(file).split(' ')
        self.start_pos = int(comment_list[1])
        self.errors_cnt = int(comment_list[2])
        self.quality = read_nonempty_line(file)


def read_nonempty_line(file):
    while True:
        raw_line = file.readline()
        if len(raw_line) == 0:
            raise EOFError()
        line = raw_line.rstrip('\n')
        if len(line) > 0:
            return line
