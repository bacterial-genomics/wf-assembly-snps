process PARSNP {
    label 'process_high'

    container "https://depot.galaxyproject.org/singularity/parsnp:2.1.3--h077b44d_0"

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

      echo "\$out" >> inputs.tmp
    done

    reference="staging/\$(ls -1S staging | head -1)"

    grep -vxF "\$reference" inputs.tmp > inputs.txt

    parsnp \
      --sequences inputs.txt \
      --validate-input \
      --reference \${reference} \
      --output-dir ./output \
      --no-maf \
      --threads $task.cpus \
      -P $mem \
      --skip-phylogeny \
      --verbose \
      $args

    harvesttools -i ./output/parsnp.ggr -M core.aln

    if [[ \$(grep -c "\\.ref" "core.aln") -eq 1 ]]; then
      echo "Stripping .ref from core genome alignment FastA"
      sed -i 's/\\.ref//1' "core.aln"
    else
      echo "'.ref' occurs multiple times or not at all; no changes made to core genome alignment FastA"
    fi
    """
}
