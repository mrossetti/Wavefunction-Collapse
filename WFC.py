'''
Credits:
    - wonderful post by Robert Heaton:
    https://robertheaton.com/2018/12/17/wavefunction-collapse-algorithm/
    - original source code here:
    https://github.com/robert/wavefunction-collapse/blob/093d218e79f01cbb7787b898b458cc5d174ad7cb/main.py#L278
'''

import colorama
from math import log
from random import random



class WFC:
    compass = {'NORTH': ( 0, -1),
               'EAST' : ( 1,  0),
               'SOUTH': ( 0,  1),
               'WEST' : (-1,  0)}

    def __init__(self, input_matrix):
        self.rules, self.freqs = self.parse(input_matrix)
        self.states = set(self.freqs.keys())

    @staticmethod
    def valid_dirs(matrix, x, y):
        hm, wm = len(matrix)-1, len(matrix[0])-1  # max index available
        dirs = []
        if y > 0: dirs.append(WFC.compass['NORTH'])
        if x <wm: dirs.append(WFC.compass['EAST'])
        if y <hm: dirs.append(WFC.compass['SOUTH'])
        if x > 0: dirs.append(WFC.compass['WEST'])
        return dirs

    def generate(self, shape, colors=None):
        states, rules, freqs = self.states, self.rules, self.freqs
        matrix = self.get_full_matrix(shape, fn_xy=lambda x,y: states.copy())
        while not self.is_all_collapsed(matrix):
            x, y = self.xy_min_entropy(matrix, freqs)
            self.collapse(matrix, x, y, freqs)
            self.propagate(matrix, x, y, rules)
        output = self.get_all_collapsed(matrix)
        if colors is not None: self.render(output, colors)
        return output

    @staticmethod
    def get_full_matrix(shape, fill_val=None, fn_xy=None):
        if fn_xy is None:
            fn_xy = lambda x,y: fill_val
        rows, cols = shape
        matrix = []
        for y in range(rows):
            row = []
            for x in range(cols):
                val = fn_xy(x, y)
                row.append(val)
            matrix.append(row)
        return matrix

    @staticmethod
    def is_all_collapsed(matrix):
        h, w = len(matrix), len(matrix[0])
        for y in range(h):
            for x in range(w):
                if len(matrix[y][x]) > 1:
                    return False
        return True

    @staticmethod
    def xy_min_entropy(matrix, freqs):
        h, w = len(matrix), len(matrix[0])
        p = freqs; noise = lambda: -random()/1000
        argmin_entropy, min_entropy = None, float('inf')
        for y in range(h):
            for x in range(w):
                possible_states = matrix[y][x]
                if len(possible_states) > 1:  # otherwise already collpased!
                    entropy = -sum([p[s]*log(p[s]) for s in possible_states])
                    noisy_entropy = entropy + noise()
                    if noisy_entropy < min_entropy:
                        min_entropy = noisy_entropy
                        argmin_entropy = (x, y)
        return argmin_entropy

    @staticmethod
    def collapse(matrix, x, y, freqs):
        # freqs contains ALL (state, freq) pairs
        # but we must pick a state only from the POSSIBLE states at this x, y
        valid_freqs = {}
        possible_states = matrix[y][x]
        total = 0
        for s in possible_states:
            valid_freqs[s] = freqs[s]
            total += freqs[s]
        # if possible states is a strict subset of all states (i.e. pos ⊂ all)
        # then valid_freqs is not normalized (i.e. Σp ≠ 1)
        chance = random() * total
        cdf = 0
        for state, prob in valid_freqs.items():
            cdf += prob
            if chance <= cdf:
                matrix[y][x] = set([state])
                return

    @staticmethod
    def propagate(matrix, x, y, rules):
        stack = [(x, y)]
        while stack:
            cur_x, cur_y = stack.pop()
            cur_possible_states = matrix[cur_y][cur_x]
            for dx, dy in WFC.valid_dirs(matrix, cur_x, cur_y):
                adj_x, adj_y = (cur_x + dx, cur_y + dy)
                adj_possible_states = matrix[adj_y][adj_x].copy()
                for adj_s in adj_possible_states:
                    adj_s_compatible = any([(cur_s, adj_s, (dx, dy)) in rules
                                            for cur_s in cur_possible_states])
                    if not adj_s_compatible:
                        matrix[adj_y][adj_x].remove(adj_s)
                        stack.append((adj_x, adj_y))

    @staticmethod
    def get_all_collapsed(matrix):
        h, w = len(matrix), len(matrix[0])
        collapsed_matrix = WFC.get_full_matrix((h, w))
        for y in range(h):
            for x in range(w):
                (state,) = matrix[y][x]
                collapsed_matrix[y][x] = state
        return collapsed_matrix

    @staticmethod
    def render(collapsed_matrix, colors):
        h, w = len(collapsed_matrix), len(collapsed_matrix[0])
        for y in range(h):
            for x in range(w):
                sym = collapsed_matrix[y][x]
                print(colors[sym] + sym + colorama.Style.RESET_ALL, end=' ')
            print()

    @staticmethod
    def parse(input_matrix):
        h, w = len(input_matrix), len(input_matrix[0])
        rules, count = set(), {}
        for y in range(h):
            for x in range(w):
                cur_tile = input_matrix[y][x]
                count[cur_tile] = count.setdefault(cur_tile, 0) + 1
                for dx, dy in WFC.valid_dirs(input_matrix, x, y):
                    adj_tile = input_matrix[y+dy][x+dx]
                    rules.add((cur_tile, adj_tile, (dx, dy)))
        tot = h*w
        freqs = {tile_sym: cnt/tot for tile_sym, cnt in count.items()}
        return rules, freqs



def main():
    input_matrix = [
        ['G','G','G','G'],
        ['G','G','G','G'],
        ['G','G','G','G'],
        ['G','S','S','G'],
        ['S','W','W','S'],
        ['W','W','W','W'],
        ['W','W','W','W'],
    ]
    colors = {
        'G': colorama.Fore.GREEN,   # Grass
        'W': colorama.Fore.CYAN,    # Water
        'S': colorama.Fore.YELLOW,  # Sand
    }
    shape = (20, 30)  # output shape (rows, columns)
    # Wavefunction Collapse
    wfc = WFC(input_matrix)
    wfc.generate(shape, colors)
    # print(wfc.freqs) {state: freq}
    # print(wfc.rules) {(state1, state2, dir)}

    # UNDERSTANDING RULES
    # The above rule is read as: f'Can {state2} stay at {dir} from {state1}'
    # To illustrate, extracted rules from the sample above (input_matrix) are:
    # - Water can stay next to water from all four dirs
    # - Sand can stay next to sand from {'EAST', 'WEST'}
    # - Grass can stay next to grass from all four dirs
    # (Water can never stay next to Grass or vice-versa in any direction)
    # - Water can stay next to sand (only) from all dirs except {'NORTH'}
    assert ('S', 'W', WFC.compass['EAST']) in wfc.rules
    assert ('S', 'W', WFC.compass['SOUTH']) in wfc.rules
    assert ('S', 'W', WFC.compass['WEST']) in wfc.rules
    assert ('S', 'W', WFC.compass['NORTH']) not in wfc.rules
    # - Sand can, therefore, stay next to water from all dirs except {'SOUTH'}
    assert ('W', 'S', WFC.compass['NORTH']) in wfc.rules
    assert ('W', 'S', WFC.compass['EAST']) in wfc.rules
    assert ('W', 'S', WFC.compass['WEST']) in wfc.rules
    assert ('W', 'S', WFC.compass['SOUTH']) not in wfc.rules
    # - Grass can stay next to sand from all dirs except {'SOUTH'}
    assert ('S', 'G', WFC.compass['NORTH']) in wfc.rules
    assert ('S', 'G', WFC.compass['EAST']) in wfc.rules
    assert ('S', 'G', WFC.compass['WEST']) in wfc.rules
    assert ('S', 'G', WFC.compass['SOUTH']) not in wfc.rules
    # - Sand can, therefore, stay next to grass from all dirs except {'NORTH'}
    assert ('G', 'S', WFC.compass['EAST']) in wfc.rules
    assert ('G', 'S', WFC.compass['SOUTH']) in wfc.rules
    assert ('G', 'S', WFC.compass['WEST']) in wfc.rules
    assert ('G', 'S', WFC.compass['NORTH']) not in wfc.rules



if __name__ == '__main__':
    main()
