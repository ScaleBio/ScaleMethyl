import nextflow.Nextflow

// Adapted from nf-core/rnaseq
class ParamLogger {
    public static void initialise(workflow, params, log) {
        if (params.help) {
            log.info("\nScaleBio Seq Suite: DNA Methylation Workflow\nFor usage information, please see README.md\n")
            System.exit(0)
        }
        log.info paramsSummaryLog(workflow, params)
        validateParams(params, log)
    }

    public static String throwError(errMessage) {
        Map colors = logColours()
        Nextflow.error("${colors.red}ERROR${colors.reset}: $errMessage")
    }
    
    // Takes a camelCase string and converts it to kebab-case
    // Necessary for validating parameters since nextflow converts camelCase to kebab-case
    public static String camelToKebab(String camelCase) {
        return camelCase.replaceAll(/([a-z])([A-Z])/, '$1-$2').toLowerCase()
    }

    // Check required workflow inputs
    public static validateParams(params, log) {
        // PLEASE ADD NEW PARAMETERS TO THE allowedParameters LIST
        def allowedParameters = [ 'runFolder', 'fastqSamplesheet', 'fastqDir', 'samples', 'libStructure', 'splitFastq',
                                  'bclConvertParams', 'genome', 'outDir', 'resultDir', 'reporting', 'fastqOut',
                                  'trimOut', 'bamOut', 'bamDedupOut', 'runTssEnrich', 'matrixGenerationCH',
                                  'matrix-generation-CH', 'matrixGenerationCG', 'matrix-generation-CG', 'bam1Dir', 'bam2Dir',
                                  'bamMergeOut', 'bamRgHeader', 'trimFastq', 'adapters', 'fastqc', 'topCellPercentile',
                                  'minCellRatio', 'minUniqCount', 'minUniqTotal', 'maxUniqTotal', 'maxMemory', 'maxCpus',
                                  'maxTime', 'dedupKey', 'minMapq', 'help', 'allcOut', 'covOut', 'chReadsThreshold', 'amethystOut' ]
        def masterListOfParams = allowedParameters
        allowedParameters.each { str ->
            masterListOfParams += camelToKebab(str)}
        def parameterDiff = params.keySet() - masterListOfParams
        if (parameterDiff.size() == 1){
            log.warn("[Argument Error] Parameter $parameterDiff is not valid in the pipeline!")
        }
        if (parameterDiff.size() > 1){
            log.warn("[Argument Error] Parameters $parameterDiff are not valid in the pipeline!")
        }
        if (params.samples == null || params.samples == true) {
            throwError("Must specify --samples (e.g. samples.csv)")
        }
        if (params.genome == null || params.genome == true) {
            throwError("Must specify --genome")
        }
        if ((params.runFolder || params.fastqDir) && (params.bam1Dir || params.bam2Dir)) {
            throwError("Either start from bam file or a runFolder/fastqDir, not both")
        }
        if (params.reporting && !params.resultDir) {
            throwError("Must specify --resultDir when --reporting is enabled")
        }
    }
    
    static LinkedHashMap paramsSummaryMap(workflow, params) {
        def Map nextflowOpts = [:] // Core Nextflow options
        def Map workflowOpts = [:] // ScaleBio Workflow parameters
        def Map execOpts = [:] // ScaleBio Workflow execution options
        def Map inputOpts = [:] // Workflow inputs

        nextflowOpts['Workflow Directory:']   = workflow.projectDir
        nextflowOpts['Workflow Version:'] = workflow.manifest.version
        if (workflow.revision) { nextflowOpts['Workflow Revision'] = workflow.revision }
        nextflowOpts['Command Line'] = workflow.commandLine
        nextflowOpts['Nextflow RunName'] = workflow.runName
        nextflowOpts['Profile'] = workflow.profile
        nextflowOpts['Config Files'] = workflow.configFiles.join(', ')
        if (workflow.containerEngine) {
            nextflowOpts['Container Engine'] = workflow.containerEngine
        }
        nextflowOpts['Launch Directory'] = workflow.launchDir
        nextflowOpts['Work Directory'] = workflow.workDir

        execOpts['Workflow Output'] = params.outDir
        if (params.reporting) { execOpts['reporting'] = params.reporting }
        execOpts['maxMemory'] = params.maxMemory
        execOpts['maxCpus'] = params.maxCpus
        execOpts['splitFastq'] = params.splitFastq
        if (params.fastqSamplesheet) { execOpts['fastqSamplesheet'] = params.fastqSamplesheet }
        execOpts['runTssEnrich'] = params.runTssEnrich
        execOpts['matrixGenerationCH'] = params.matrixGenerationCH
        execOpts['matrixGenerationCG'] = params.matrixGenerationCG
        execOpts['trimFastq'] = params.trimFastq

        if (params.fastqDir) {
            inputOpts['fastqDir'] = params.fastqDir
        } else if (params.runFolder) {
            inputOpts['runFolder'] = params.runFolder
        } else if (params.bam1Dir && params.bam2Dir) {
            inputOpts['bam1Dir'] = params.bam1Dir
            inputOpts['bam2Dir'] = params.bam2Dir
        }
        if (params.reporting) { inputOpts['resultDir'] = params.resultDir }
        inputOpts['samples'] = params.samples
        inputOpts['genome'] = params.genome
        inputOpts['libStructure'] = params.libStructure

        workflowOpts['adapters'] = params.adapters
        workflowOpts['dedupKey'] = params.dedupKey
        workflowOpts['minMapq'] = params.minMapq
        workflowOpts['bclConvertParams'] = params.bclConvertParams
        workflowOpts['minUniqTotal'] = params.minUniqTotal
        workflowOpts['maxUniqTotal'] = params.maxUniqTotal
        workflowOpts['topCellPercentile'] = params.topCellPercentile

        return [ 'Core Nextflow Options': nextflowOpts,
                 'Inputs': inputOpts,
                 'Workflow Execution Options': execOpts,
                 'Analysis Parameters': workflowOpts
        ]
    }

    // Beautify parameters for summary and return as string
    public static String paramsSummaryLog(workflow, params) {
        Map colors = logColours()
        String output  = ''
        def paramsMap = paramsSummaryMap(workflow, params)
        def maxChars  = paramsMaxChars(paramsMap)
        for (group in paramsMap.keySet()) {
            def groupParams = paramsMap.get(group)  // This gets the parameters of that particular group
            if (groupParams) {
                output += colors.bold + group + colors.reset + '\n'
                for (param in groupParams.keySet()) {
                    output += "  " + colors.blue + param.padRight(maxChars) + ": " + colors.green +  groupParams.get(param) + colors.reset + '\n'
                }
                output += '\n'
            }
        }
        output += dashedLine()
        return output
    }

    public static Map logColours(monochromeLogs=false) {
        Map colorcodes = [:]
        // Reset / Meta
        colorcodes['reset']      = monochromeLogs ? '' : "\033[0m"
        colorcodes['bold']       = monochromeLogs ? '' : "\033[1m"
        colorcodes['dim']        = monochromeLogs ? '' : "\033[2m"
        colorcodes['underlined'] = monochromeLogs ? '' : "\033[4m"
        colorcodes['blink']      = monochromeLogs ? '' : "\033[5m"
        colorcodes['reverse']    = monochromeLogs ? '' : "\033[7m"
        colorcodes['hidden']     = monochromeLogs ? '' : "\033[8m"

        // Regular Colors
        colorcodes['black']      = monochromeLogs ? '' : "\033[0;30m"
        colorcodes['red']        = monochromeLogs ? '' : "\033[0;31m"
        colorcodes['green']      = monochromeLogs ? '' : "\033[0;32m"
        colorcodes['yellow']     = monochromeLogs ? '' : "\033[0;33m"
        colorcodes['blue']       = monochromeLogs ? '' : "\033[0;34m"
        colorcodes['purple']     = monochromeLogs ? '' : "\033[0;35m"
        colorcodes['cyan']       = monochromeLogs ? '' : "\033[0;36m"
        colorcodes['white']      = monochromeLogs ? '' : "\033[0;37m"

        // Bold
        colorcodes['bblack']     = monochromeLogs ? '' : "\033[1;30m"
        colorcodes['bred']       = monochromeLogs ? '' : "\033[1;31m"
        colorcodes['bgreen']     = monochromeLogs ? '' : "\033[1;32m"
        colorcodes['byellow']    = monochromeLogs ? '' : "\033[1;33m"
        colorcodes['bblue']      = monochromeLogs ? '' : "\033[1;34m"
        colorcodes['bpurple']    = monochromeLogs ? '' : "\033[1;35m"
        colorcodes['bcyan']      = monochromeLogs ? '' : "\033[1;36m"
        colorcodes['bwhite']     = monochromeLogs ? '' : "\033[1;37m"

        // Underline
        colorcodes['ublack']     = monochromeLogs ? '' : "\033[4;30m"
        colorcodes['ured']       = monochromeLogs ? '' : "\033[4;31m"
        colorcodes['ugreen']     = monochromeLogs ? '' : "\033[4;32m"
        colorcodes['uyellow']    = monochromeLogs ? '' : "\033[4;33m"
        colorcodes['ublue']      = monochromeLogs ? '' : "\033[4;34m"
        colorcodes['upurple']    = monochromeLogs ? '' : "\033[4;35m"
        colorcodes['ucyan']      = monochromeLogs ? '' : "\033[4;36m"
        colorcodes['uwhite']     = monochromeLogs ? '' : "\033[4;37m"

        // High Intensity
        colorcodes['iblack']     = monochromeLogs ? '' : "\033[0;90m"
        colorcodes['ired']       = monochromeLogs ? '' : "\033[0;91m"
        colorcodes['igreen']     = monochromeLogs ? '' : "\033[0;92m"
        colorcodes['iyellow']    = monochromeLogs ? '' : "\033[0;93m"
        colorcodes['iblue']      = monochromeLogs ? '' : "\033[0;94m"
        colorcodes['ipurple']    = monochromeLogs ? '' : "\033[0;95m"
        colorcodes['icyan']      = monochromeLogs ? '' : "\033[0;96m"
        colorcodes['iwhite']     = monochromeLogs ? '' : "\033[0;97m"

        // Bold High Intensity
        colorcodes['biblack']    = monochromeLogs ? '' : "\033[1;90m"
        colorcodes['bired']      = monochromeLogs ? '' : "\033[1;91m"
        colorcodes['bigreen']    = monochromeLogs ? '' : "\033[1;92m"
        colorcodes['biyellow']   = monochromeLogs ? '' : "\033[1;93m"
        colorcodes['biblue']     = monochromeLogs ? '' : "\033[1;94m"
        colorcodes['bipurple']   = monochromeLogs ? '' : "\033[1;95m"
        colorcodes['bicyan']     = monochromeLogs ? '' : "\033[1;96m"
        colorcodes['biwhite']    = monochromeLogs ? '' : "\033[1;97m"

        return colorcodes
    }

    // Creates a dashed line 
    public static String dashedLine() {
        Map colors = logColours()
        return "-${colors.dim}----------------------------------------------------${colors.reset}-"
    }
    private static Integer paramsMaxChars(paramsMap) {
        Integer maxChars = 0
        for (group in paramsMap.keySet()) {
            def groupParams = paramsMap.get(group)  // This gets the parameters of that particular group
            for (param in groupParams.keySet()) {
                if (param.size() > maxChars) {
                    maxChars = param.size()
                }
            }
        }
        return maxChars
    }
}
