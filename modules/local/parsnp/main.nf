process PARSNP {
    label 'process_high'

    container "https://depot.galaxyproject.org/singularity/parsnp:2.1.5--h077b44d_0"

    input:
    path(fasta, stageAs: "input/*")
    path(reference)

    output:
    path("core.aln"), emit: aln

    when:
    task.ext.when == null || task.ext.when

    script:
    def args      = task.ext.args   ?: ''
    //prefix        = task.ext.prefix ?: "${meta.id}"
    def mem = ((task.memory as MemoryUnit).toBytes() * 0.8) as int

    """
    mkdir staging

    shopt -s nullglob

    for f in input/*; do
      [[ -f "\$f" ]] || continue

      base="\$(basename "\$f")"

      if [[ "\$f" == *.gz ]]; then
        out="staging/\${base%.gz}"
        gunzip -c "\$f" > "\$out"
      else
        out="staging/\$base"
        cp "\$f" "\$out"
      fi
    done

    reference="\$(ls -1S staging | head -1)"

    mv staging/\$reference ./\$reference

    parsnp \
      --sequences ./staging \
      --validate-input \
      --reference ./\$(basename \$reference) \
      --output-dir ./output \
      --curated \
      --no-maf \
      --threads $task.cpus \
      -P $mem \
      --skip-phylogeny \
      --verbose \
      $args

    harvesttools -x output/parsnp.xmfa -M core.aln

    if [[ \$(grep -c "\\.ref" "core.aln") -eq 1 ]]; then
      echo "Stripping .ref from core genome alignment FastA"
      sed -i 's|\\.ref||g' ./core.aln
    else
      echo "'.ref' occurs multiple times or not at all; no changes made to core genome alignment FastA"
    fi
    """
}
