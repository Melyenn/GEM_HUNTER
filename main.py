import time
import copy
import os
import sys
from itertools import combinations, product
from pysat.solvers import Glucose3


def read_file(file_path):
    with open(file_path, "r", encoding='utf-8') as f:
        return [line.strip().split(", ") for line in f.readlines()]


def get_neighbors(grid, x, y, row, col):
    neighbors = []
    dx = [0, 0, -1, 1, -1, -1, 1, 1]
    dy = [-1, 1, 0, 0, -1, 1, -1, 1]

    for i in range(8):
        nx, ny = x + dx[i], y + dy[i]
        if 0 <= nx < row and 0 <= ny < col and not grid[nx][ny].isdigit():
            neighbors.append((nx, ny))
    return neighbors


def generate_cnf(grid, row, col):
    clauses = set()
    for x in range(row):
        for y in range(col):
            if grid[x][y].isdigit():
                num_traps = int(grid[x][y])

                neighbors_trap = get_neighbors(grid, x, y, row, col)
                num_emptys = len(neighbors_trap)

                if num_traps < 1 or num_traps > num_emptys:
                    return None

                if num_emptys >= num_traps:
                    for subset in combinations(neighbors_trap, num_emptys - num_traps + 1):
                        clause = tuple(sorted(nx * col + ny + 1 for (nx, ny) in subset))
                        clauses.add(clause)

                if num_emptys > num_traps:
                    for subset in combinations(neighbors_trap, num_traps + 1):
                        clause = tuple(sorted(-(nx * col + ny + 1) for (nx, ny) in subset))
                        clauses.add(clause)

    return [list(clause) for clause in clauses]


def solve_with_pysat(grid, empty_cells, col, clauses):
    solver = Glucose3()

    for clause in clauses:
        solver.add_clause(clause)

    if solver.solve():
        answer = solver.get_model()
        if answer is None:
            return None

        result = [row.copy() for row in grid]

        for x, y in empty_cells:
            result[x][y] = 'T' if x * col + y + 1 in answer else 'G'
        return result

    return None


def write_file(grid, filepath):
    with open(filepath, "w", encoding='utf-8') as f:
        for row in grid:
            print(', '.join(str(elem) for elem in row), file=f)


def check_condition_stop(grid, row, col, cnf_clauses, empty_cells):
    ans = ['_'] * (row * col + 1)
    for x, y in empty_cells:
        ans[x * col + y + 1] = grid[x][y]

    for clause in cnf_clauses:
        check = False
        for c in clause:
            if (c < 0 and ans[abs(c)] == 'G') or (c > 0 and ans[c] == 'T'):
                check = True
                break
        if not check:
            return False
    return True


def solve_bruteforce(grid, row, col, cnf_clauses, empty_cells):
    test_grid = [row.copy() for row in grid]
    for assignment in product(['G', 'T'], repeat=len(empty_cells)):
        for (x, y), val in zip(empty_cells, assignment):
            test_grid[x][y] = val
        if check_condition_stop(test_grid, row, col, cnf_clauses, empty_cells):
            return test_grid
    return None

def unit_propagation(grid, col, cnf_clauses):
    isvalid = set()
    i = 0
    while i < len(cnf_clauses):
        clause = cnf_clauses[i]
        if len(clause) == 1:
            c = clause[0]
            if -c in isvalid:
                return None
            if c > 0:
                x, y = (c - 1) // col, (c - 1) % col
                grid[x][y] = 'T'
                isvalid.add(c)
            else:
                c = abs(c)
                x, y = (c - 1) // col, (c - 1) % col
                grid[x][y] = 'G'
                isvalid.add(-c)
            cnf_clauses.remove(clause)
        else:
            i += 1
    for val in isvalid:
        i = 0
        while i < len(cnf_clauses):
            if -val in cnf_clauses[i]:
                cnf_clauses[i].remove(-val)
            if val in cnf_clauses[i]:
                cnf_clauses.remove(cnf_clauses[i])
            else:
                i += 1

    return [list(clause) for clause in set(tuple(sorted(clause)) for clause in cnf_clauses)]

# def pure_literal_elimination(grid, col, cnf_clauses):
#     if cnf_clauses is None:
#         return None
#     abs_val = [abs(val) for clause in cnf_clauses for val in clause]
#     all_val = [val for clause in cnf_clauses for val in clause]
#     all_valid = []
#     for val in abs_val:
#         if val in all_val and -val not in all_val:
#             x, y = (val - 1) // col, (val - 1) % col
#             grid[x][y] = 'T'
#             all_valid.append(val)
#         if -val in all_val and val not in all_val:
#             x, y = (val - 1) // col, (val - 1) % col
#             grid[x][y] = 'G'
#             all_valid.append(-val)
#     for val in all_valid:
#         i = 1
#         while i < len(cnf_clauses):
#             if val in cnf_clauses[i]:
#                 cnf_clauses.remove(cnf_clauses)
#             else:
#                 i += 1
#     return cnf_clauses

def solve_dpll(grid, col, idx, empty_cells, cnf_clauses):
    cnf_clauses = unit_propagation(grid, col, cnf_clauses)
    # cnf_clauses = pure_literal_elimination(grid, col, cnf_clauses)

    if cnf_clauses is None:
        return None

    if len(cnf_clauses) == 0:
        while idx < len(empty_cells):
            x, y = empty_cells[idx]
            if grid[x][y] == '_':
                grid[x][y] = 'G'
            idx += 1
        return grid

    if any(len(clause) == 0 for clause in cnf_clauses):
        return None

    while idx < len(empty_cells):
        x, y = empty_cells[idx]
        if grid[x][y] == '_':
            val = x * col + y + 1
            break
        idx += 1
    else:
        return None

    new_cnf = copy.deepcopy(cnf_clauses)
    new_cnf.append([val])
    new_grid = [row.copy() for row in grid]
    new_grid = solve_dpll(new_grid, col, idx + 1, empty_cells, new_cnf)
    if new_grid:
        return new_grid

    new_cnf = copy.deepcopy(cnf_clauses)
    new_cnf.append([-val])
    new_grid = [row.copy() for row in grid]
    return solve_dpll(new_grid, col, idx + 1, empty_cells, new_cnf)


def main():
    if len(sys.argv) < 4:
        print('Usage: python main.py <num_method> <method> <input_file>')
        sys.exit(1)
    num_method = int(sys.argv[1])
    idx = 2
    methods = []
    while idx < num_method + 2:
        if sys.argv[idx] not in {'--pysat', '--bruteforce', '--backtracking'}:
            print('Error: Invalid method')
            sys.exit(1)
        methods.append(sys.argv[idx])
        idx += 1

    input_file = sys.argv[idx]
    grid = read_file(input_file)
    rows, cols = len(grid), len(grid[0])
    empty_cells = [(x, y) for x in range(rows) for y in range(cols) if grid[x][y] == '_']

    start_time = time.time()
    cnf_clauses = generate_cnf(grid, rows, cols)
    end_time = time.time()
    print(f'Generate CNF time: {end_time - start_time}s')
    if not cnf_clauses:
        print('No valid CNF found')
        sys.exit(1)

    convert_methods = {
        '--pysat': solve_with_pysat,
        '--bruteforce': solve_bruteforce,
        '--backtracking': solve_dpll
    }

    for method in methods:
        solver = convert_methods[method]
        os.makedirs(f'testcases/{method[2:]}', exist_ok=True)
        output_file = input_file.replace('testcases/', f'testcases/{method[2:]}/')
        output_file = output_file.replace('input_', 'output_')
        print(f'Solving with {method[2:]}...')
        start_time = time.time()
        match method:
            case '--backtracking':
                result = solver(grid, cols, 0, empty_cells, copy.deepcopy(cnf_clauses))
            case '--bruteforce':
                result = solver(grid, rows, cols, cnf_clauses, empty_cells)
            case '--pysat':
                result = solver(grid, empty_cells, cols, cnf_clauses)
            case _:
                print("Error!")
                exit(0)
        print(f'{method[2:]} time: {time.time() - start_time: }s')
        if result:
            print(f'Write result in file: {output_file}')
            write_file(result, output_file)
        else:
            print('No solution')

if __name__ == '__main__':
    main()
