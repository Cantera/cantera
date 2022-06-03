import textwrap

i = 0

def normalize(code: str, indent: int = 0, trim_first: bool = False):
    code = textwrap.dedent(code).strip()

    if indent:
        code = textwrap.indent(code, ' ' * indent)
    
        if trim_first:
            code = code[indent:]

    return code
