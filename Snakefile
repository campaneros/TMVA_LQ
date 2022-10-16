rule copy:
    input: 'input.txt'
    output: 'output.txt'
    shell: 'cp {input} {output}'
