class: CommandLineTool
cwlVersion: v1.0
$namespaces:
  sbg: 'https://www.sevenbridges.com/'
id: samtools_coverage
baseCommand:
  - samtools
  - coverage
inputs:
  - id: include_secondary
    type: boolean?
    default: false
    doc: 'If true, include SECONDARY alignments by excluding only UNMAP,QCFAIL,DUP (omit SECONDARY from default exclusions).'
  - id: input
    type: File
    inputBinding:
      position: 1
arguments:
  - prefix: --excl-flags
    valueFrom: '$(inputs.include_secondary ? "UNMAP,QCFAIL,DUP" : "UNMAP,SECONDARY,QCFAIL,DUP")'
outputs:
  - id: output
    type: File
    outputBinding:
      glob: virus.coverage.txt
label: samtools-coverage
requirements:
  - class: DockerRequirement
    dockerPull: 'quay.io/biocontainers/samtools:1.15--h1170115_1'
  - class: InlineJavascriptRequirement
stdout: virus.coverage.txt
