"""
Credits
    - read this wonderful post by Robert Heaton:
    https://robertheaton.com/2018/12/17/wavefunction-collapse-algorithm/
    - find the original source code here:
    https://github.com/robert/wavefunction-collapse/blob/093d218e79f01cbb7787b898b458cc5d174ad7cb/main.py#L278
"""

from math import log
from random import random, choices
from collections import defaultdict
from itertools import product
import colorama


def render(collapsed_matrix, states_colors={}, states_symbols={}):
    colors = list(vars(colorama.Fore).values())
    h, w = len(collapsed_matrix), len(collapsed_matrix[0])

    for y in range(h):
        for x in range(w):
            s = collapsed_matrix[y][x]
            color = states_colors.setdefault(s, choices(colors)[0])
            symbol = states_symbols.get(s, s)
            print(color + symbol + colorama.Style.RESET_ALL, end=' ')
        print()


def directions_valid(matrix, x, y):
    """
    return all (dx, dy) available (in-bounds) from x, y
    """
    row_max_idx = len(matrix)-1
    col_max_idx = len(matrix[0])-1
    
    directions = []
    
    if x > 0:
        directions.append((-1, 0))  # west
    if y > 0:
        directions.append((0, -1))  # north
    if x < col_max_idx:
        directions.append((1, 0))   # east
    if y < row_max_idx:
        directions.append((0, 1))   # south
    
    return directions


class WavefunctionCollapse:

    def extract(self, input_matrix):
        """
        rules: {(tile state, adjacent tile state, direction from tile observed), ...}
        freqs: {tile state: frequency}
        return rules, freqs
        """
        rules, count = set(), defaultdict(int)

        h, w = len(input_matrix), len(input_matrix[0])
        
        for x, y in product(range(w), range(h)):
            cur_tile = input_matrix[y][x]
            count[cur_tile] += 1
            
            for dx, dy in directions_valid(input_matrix, x, y):
                adj_tile = input_matrix[y + dy][x + dx]
                rules.add((cur_tile, adj_tile, (dx, dy)))
        
        total = h * w
        freqs = {state: count / total for state, count in count.items()}
        
        self._input_matrix = input_matrix
        self.rules = rules
        self.freqs = freqs
    
    def generate(self, output_matrix_shape=None):
        """
        Initially, any tile is potentially in all states.
        (every tile has the same entropy)
        
        WFC iteration:
            - Collpase() tile state (noisy best first search min entropy)
            - Propagate() collapsed neighbours (using rules)
        
        The tile with minimum entropy is that tile whose #possible_states
        (according to currently collpased neighbors and rules) is minimum.
        
        WFC may reach a point in which a tile cannot collapse (i.e. zero
        possible states). This paradox forces WFC to completely restart
        (a better implementation could make use of backtracking).
        """
        h, w = output_matrix_shape or self._input_matrix.shape
        possible_states = lambda: set(self.freqs.keys())
        matrix = [[possible_states() for x in range(w)] for y in range(h)]
        
        while not self.is_all_collapsed(matrix):
            x, y = self.xy_min_entropy(matrix, self.freqs)
            self.collapse(matrix, x, y, self.freqs)
            self.propagate(matrix, x, y, self.rules)
        
        output = self.get_all_collapsed(matrix)
        
        return output
    
    @staticmethod
    def is_all_collapsed(matrix):
        return all(len(possible_states) == 1
                   for row in matrix for possible_states in row)
    
    @staticmethod
    def get_all_collapsed(matrix):
        return [[S.pop() for S in row] for row in matrix]

    @staticmethod
    def xy_min_entropy(matrix, p):
        """
        return x, y of tile with lowest entropy (with the min #possible_states > 1)
        """
        # possible_states in a given x, y (S)
        h, w = len(matrix), len(matrix[0])
        H = lambda S: -sum([p[s] * log(p[s]) for s in S])
        noise = lambda: -random() / 1000.0  # settle ties

        minH, xy = min( ((H(matrix[y][x]) + noise(), (x, y))
                        for x, y in product(range(w), range(h))
                        if len(matrix[y][x]) > 1) )
        
        return xy

    @staticmethod
    def collapse(matrix, x, y, freqs):
        """
        modifies matrix in-place at (x, y) index collpasing the tile to a single
        state (among the possible ones). Each (possible) state's probability of
        being selected is proportional to its frequency.
        """
        possible_states = list(matrix[y][x])
        probs = [freqs[s] for s in possible_states]
        matrix[y][x] = choices(possible_states, probs)

    @staticmethod
    def propagate(matrix, x, y, rules):
        """
        modifies matrix in-place by enforcing rules starting from collapsed x, y.
        Is there any adj_s not compatible that we should remove from the possible
        states in this (adj_y, adj_x) location?
        if a specific adj_s (among the possible states in adj) is not compatible
        (under cur-adj orientation) with any of the possible states in cur (cur_s)
        then, there is no way we could every use adj_s in adj, hence we remove it
        from the possible states in adj_s. This partial/full collapse in adj_s can
        however result in other rules propagation, hence propagate will recursively
        check its adjacents as well. 
        """
        stack = [(x, y)]
        
        while stack:
            cur_x, cur_y = stack.pop()
            cur_possible_states = matrix[cur_y][cur_x]
            
            for dx, dy in directions_valid(matrix, cur_x, cur_y):
                adj_x, adj_y = (cur_x + dx, cur_y + dy)
                adj_possible_states = matrix[adj_y][adj_x].copy()

                for adj_s in adj_possible_states:
                    adj_s_compatible = any([(cur_s, adj_s, (dx, dy)) in rules
                                            for cur_s in cur_possible_states])
                    
                    if not adj_s_compatible:
                        matrix[adj_y][adj_x].remove(adj_s)
                
                # if incompatibility found
                if len(matrix[adj_y][adj_x]) < len(adj_possible_states):
                    stack.append((adj_x, adj_y))






def main():
    # Input from which to extract rules
    input_matrix = [
        ['G','G','G','G'],
        ['G','G','G','G'],
        ['G','G','G','G'],
        ['G','S','S','G'],
        ['S','W','W','S'],
        ['W','W','W','W'],
        ['W','W','W','W'],
    ]
    
    # Wavefunction Collapse
    wfc = WavefunctionCollapse()
    wfc.extract(input_matrix)

    # print(wfc.freqs)  # {state: freq}
    # print(wfc.rules)  # {(state1, state2, direction)}
                        # f'Can {state2} stay at {direction} from {state1}'

    compass = {'NORTH': ( 0, -1),
               'EAST' : ( 1,  0),
               'SOUTH': ( 0,  1),
               'WEST' : (-1,  0)}
    # To illustrate, extracted rules from the sample above (input_matrix) are:
    # - Water can stay next to water from all four directions
    # - Sand can stay next to sand from {'EAST', 'WEST'}
    # - Grass can stay next to grass from all four directions
    # (Water can never stay next to Grass or vice-versa in any direction)
    # - Water can stay next to sand (only) from all directions except {'NORTH'}
    assert ('S', 'W', compass['EAST']) in wfc.rules
    assert ('S', 'W', compass['SOUTH']) in wfc.rules
    assert ('S', 'W', compass['WEST']) in wfc.rules
    assert ('S', 'W', compass['NORTH']) not in wfc.rules
    # - Sand can, therefore, stay next to water from all dirs except {'SOUTH'}
    assert ('W', 'S', compass['NORTH']) in wfc.rules
    assert ('W', 'S', compass['EAST']) in wfc.rules
    assert ('W', 'S', compass['WEST']) in wfc.rules
    assert ('W', 'S', compass['SOUTH']) not in wfc.rules
    # - Grass can stay next to sand from all dirs except {'SOUTH'}
    assert ('S', 'G', compass['NORTH']) in wfc.rules
    assert ('S', 'G', compass['EAST']) in wfc.rules
    assert ('S', 'G', compass['WEST']) in wfc.rules
    assert ('S', 'G', compass['SOUTH']) not in wfc.rules
    # - Sand can, therefore, stay next to grass from all dirs except {'NORTH'}
    assert ('G', 'S', compass['EAST']) in wfc.rules
    assert ('G', 'S', compass['SOUTH']) in wfc.rules
    assert ('G', 'S', compass['WEST']) in wfc.rules
    assert ('G', 'S', compass['NORTH']) not in wfc.rules

    shp = (20, 30)  # desired output shape (rows, columns)
    out = wfc.generate(shp)

    colors = {
        'G': colorama.Fore.GREEN,   # Grass
        'W': colorama.Fore.CYAN,    # Water
        'S': colorama.Fore.YELLOW,  # Sand
    }

    render(out, colors)



if __name__ == '__main__':
    main()
