//
// This file holds several functions specific to the workflow/blast.nf in the wf-assembly-snps pipeline
//

import groovy.text.SimpleTemplateEngine

class WorkflowSNPS {

    //
    // Check and validate parameters
    //
    public static void initialise(params, log) {
        if (!params.input) {
            log.error "Input directory is required to perform analysis."
            System.exit(1)
        }
    }
}
