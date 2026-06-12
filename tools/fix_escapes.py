r"""Fix invalid escape sequences (SyntaxWarning in Python 3.12) by doubling
the backslash inside string literals, e.g. '$\mu$' -> '$\\mu$'.
The runtime string value is identical; verified by comparing ASTs."""
import ast
import io
import re
import sys
import tokenize

# Characters that form a valid escape sequence after a backslash
VALID = set("\n\r\\'\"abfnrtv01234567xNuU")


def fix_literal(tok_str):
    m = re.match(r"([A-Za-z]*)('''|\"\"\"|'|\")", tok_str)
    if m is None or 'r' in m.group(1).lower():
        return tok_str
    body_start = m.end(1)
    out = list(tok_str[:body_start])
    i = body_start
    while i < len(tok_str):
        c = tok_str[i]
        if c == '\\' and i + 1 < len(tok_str):
            nxt = tok_str[i + 1]
            if nxt in VALID:
                out.append(c)
            else:
                out.append('\\\\')
            out.append(nxt)
            i += 2
        else:
            out.append(c)
            i += 1
    return ''.join(out)


def main(path):
    with open(path, encoding='utf-8') as f:
        src = f.read()
    lines = src.splitlines(keepends=True)

    edits = []  # (srow, scol, erow, ecol, replacement)
    for tok in tokenize.generate_tokens(io.StringIO(src).readline):
        if tok.type == tokenize.STRING:
            fixed = fix_literal(tok.string)
            if fixed != tok.string:
                edits.append((tok.start, tok.end, fixed))

    if not edits:
        print(f'{path}: nothing to fix')
        return

    for (srow, scol), (erow, ecol), fixed in reversed(edits):
        if srow == erow:
            line = lines[srow - 1]
            lines[srow - 1] = line[:scol] + fixed + line[ecol:]
        else:
            head = lines[srow - 1][:scol]
            tail = lines[erow - 1][ecol:]
            lines[srow - 1:erow] = [head + fixed + tail]

    new_src = ''.join(lines)
    if ast.dump(ast.parse(src)) != ast.dump(ast.parse(new_src)):
        raise SystemExit(f'{path}: AST changed, aborting!')
    with open(path, 'w', encoding='utf-8', newline='') as f:
        f.write(new_src)
    print(f'{path}: fixed {len(edits)} string literals (AST verified identical)')


if __name__ == '__main__':
    for p in sys.argv[1:]:
        main(p)
