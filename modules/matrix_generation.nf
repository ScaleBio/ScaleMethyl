/**
* Module that generates merged matrices with methylation calls for all barcodes
*
* Process:
*     MtxCG
*     MtxCH
*     MergeMtxCG
*     MergeMtxCH
*     CreateALLC
*     CreateBismark
*/

// Generate binned matrix for each barcode
process MtxCG {
input:
	path(tiles) // Non overlapping genome bins
	tuple val(sample), val(sample_with_well_coordinate), path(allCells), path(metCG)
output:
	tuple val(sample_with_well_coordinate), path("${sample_with_well_coordinate}.CG.score.mtx.gz"), path("${sample_with_well_coordinate}.barcodes.tsv"), path("${sample_with_well_coordinate}.CG.features.tsv"), emit: mtx, optional: true
tag "$sample_with_well_coordinate"
label 'pyProcess'
script:
"""
	create_mtx.py --bedfile $tiles --met_calls $metCG --all_cells $allCells --sample $sample_with_well_coordinate
"""
}

// Generate binned matrix for each barcode
process MtxCH {
input:
	path(tiles) // Non overlapping genome bins
	tuple val(sample), val(sample_with_well_coordinate), path(allCells), path(metCH)
output:
	tuple val(sample_with_well_coordinate), path("${sample_with_well_coordinate}.CH.mtx.gz"), path("${sample_with_well_coordinate}.barcodes.tsv"), path("${sample_with_well_coordinate}.CH.features.tsv"), emit: mtx, optional: true
tag "$sample_with_well_coordinate"
label 'pyProcess'
script:
"""
	create_mtx.py --bedfile $tiles --met_calls $metCH --all_cells $allCells --sample $sample_with_well_coordinate --CH
"""
}

// Generate merged CG matrix with all barcodes
process MergeMtxCG {
input:
	tuple val(sample), path(scoreMtx), path(barcodes), path(features)
output:
	tuple val(sample), path("${sample}.CG.score.mtx.gz")
	tuple val(sample), path("${sample}.barcodes.tsv")
	tuple val(sample), path("${sample}.CG.features.tsv")
tag "$sample"
label 'pyProcess'
publishDir file(params.outDir) / "samples" / "genome_bin_matrix", pattern: "*.{mtx.gz,tsv}", mode: 'copy'
script:
"""
	merge_mtx.py --sample ${sample} --barcodes ${barcodes} --features ${features[0]} --scoreMtx $scoreMtx
"""
}

// Generate merged CH matrix with all barcodes
process MergeMtxCH {
input:
	tuple val(sample), path(scoreMtx), path(barcodes), path(features)
output:
	tuple val(sample), path("${sample}.CH.mtx.gz")
	tuple val(sample), path("${sample}.barcodes.tsv")
	tuple val(sample), path("${sample}.CH.features.tsv")
tag "$sample"
label 'pyProcess'
publishDir file(params.outDir) / "samples" / "genome_bin_matrix", pattern: "*.{mtx.gz,tsv}", mode: 'copy'
script:
"""
	merge_mtx.py --sample ${sample} --barcodes ${barcodes} --features ${features[0]} --scoreMtx $scoreMtx
"""
}

process CreateALLC {
input:
	tuple val(sample), val(sample_with_well_coordinate), path(allCells), path(met)
output:
	tuple val(sample), path("{CG,CH}/*.allc.tsv.gz"), optional: true
tag "$sample"
label 'pyProcess'
publishDir { outDir }, pattern: "{CG,CH}/*.allc.tsv.gz", mode: 'copy'
script:
    outDir = file(params.outDir) / "samples" / "methylation_coverage" / "allc" / sample.tokenize('.')[0]
"""
	write_allc.py --met_calls $met --all_cells $allCells --sample $sample_with_well_coordinate
"""
}

process CreateBismark {
input:
	tuple val(sample), val(sample_with_well_coordinate), path(allCells), path(met)
output:
	tuple val(sample), path("{CG,CH}/*.cov.gz"), optional: true
tag "$sample"
label 'pyProcess'
publishDir { outDir }, pattern: "{CG,CH}/*.cov.gz", mode: 'copy'
script:
    outDir = file(params.outDir) / "samples" / "methylation_coverage" / "cov" / sample.tokenize('.')[0]
"""
	write_bismark.py --met_calls $met --all_cells $allCells --sample $sample_with_well_coordinate
"""
}


workflow MATRIX_GENERATION {
take:
    allCells // Information about all cell barcodes
    metCG // Per TN5 CG methylation stats
    metCH // Per TN5 CH methylation stats
    genome // Reference genome

main:
    if (params.matrixGenerationCG) {
        metCG
            .map {
                // Retrieve sample name without well coordinate to join with allCells
                def sample = it[0].tokenize('.')[0]
                return tuple(sample, it[0], it[1])
            }
            .dump(tag: 'metCG')
            .set { metCG }
        allCells.cross(metCG)
            .map {
                // sample name, sample name.well coordinate, allCells, met calls
                tuple(it[0][0], it[1][1], it[0][1], it[1][2])
            }
            .dump(tag: 'cgPass')
            .set { cgPass }
        MtxCG(genome.genomeTiles, cgPass)
        MtxCG.out.mtx
            .dump(tag: 'mtxCG')
            // combine the outputs per sample
            .map { id, scoreMtx, barcodes, features ->
                sample = id.tokenize('.')[0]
                tuple sample, scoreMtx, barcodes, features
            }
            .groupTuple()
            .map { sample, scoreMtx, barcodes, features ->
                // flatten the list of files for each well coordinate
                tuple sample, scoreMtx.flatten(), barcodes.flatten(), features.flatten()
            }
            .dump(tag: 'mtxCGPerSample')
            .set { mtxCGPerSample }
        MergeMtxCG(mtxCGPerSample)
    }
    
    if (params.matrixGenerationCH) {
        metCH
            .map {
                // Retrieve sample name without well coordinate to join with allCells
                def sample = it[0].tokenize('.')[0]
                return tuple(sample, it[0], it[1])
            }
            .dump(tag: 'metCH')
            .set { metCH }
        allCells.cross(metCH)
            .map {
                // sample name, sample name.well coordinate, allCells, met calls
                tuple(it[0][0], it[1][1], it[0][1], it[1][2])
            }
            .dump(tag: 'chPass')
            .set { chPass }
        MtxCH(genome.genomeTilesCh, chPass)
        MtxCH.out.mtx
            .dump(tag: 'mtxCH')
            // combine the outputs per sample
            .map { id, scoreMtx, barcodes, features ->
                sample = id.tokenize('.')[0]
                tuple sample, scoreMtx, barcodes, features
            }
            .groupTuple()
            .map { sample, scoreMtx, barcodes, features ->
                // flatten the list of files for each well coordinate
                tuple sample, scoreMtx.flatten(), barcodes.flatten(), features.flatten()
            }
            .dump(tag: 'mtxCHPerSample')
            .set { mtxCHPerSample }
        MergeMtxCH(mtxCHPerSample)
    }
    if (params.allcOut) {
        CreateALLC(cgPass.concat(chPass))
    }
    if (params.covOut) {
        CreateBismark(cgPass.concat(chPass))
    }
}