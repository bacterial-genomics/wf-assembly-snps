process PARSNP {
    //tag "$meta.id"
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

    """
    mkdir staging

    for f in input/*.gz; do
      base="\$(basename "\${f%.gz}")"
      out="staging/\${base%.*}"
      gunzip -c "\$f" > "\$out"
      echo "\$out" >> inputs.txt
    done

    parsnp \
      --sequences inputs.txt \
      --reference $reference \
      --output-dir ./output \
      --threads $task.cpus \
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
