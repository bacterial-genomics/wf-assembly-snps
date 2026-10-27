process PARSNP {
    label 'process_medium'

    container "quay.io/biocontainers/parsnp:2.1.6--h3be2455_0"

    input:
    path(fasta, stageAs: "input/*")
    path(reference)

    output:
    path("core.aln"), emit: aln
    path("parsnp.ggr"), emit: ggr
    path("parsnp.log"), emit: log
    path("parsnp.xmfa"), emit: xmfa

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

    reference="\$(ls -1SL staging | head -1)"

    mv -f staging/\$reference ./\$reference

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

    mv -f output/log/parsnpAligner.log parsnp.log
    mv -f output/parsnp.ggr parsnp.ggr
    mv -f output/parsnp.xmfa parsnp.xmfa

    harvesttools -x parsnp.xmfa -M core.aln

    if [[ \$(grep -c "\\.ref" "core.aln") -eq 1 ]]; then
      echo "Stripping .ref from core genome alignment FastA"
      sed -i 's|\\.ref||1' core.aln
    else
      echo "'.ref' occurs multiple times or not at all; no changes made to core genome alignment FastA"
    fi
    """
}
