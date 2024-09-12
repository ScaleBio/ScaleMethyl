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
*     CreateAmethyst
*/

// Generate binned matrix for each barcode
process MtxCG {
input:
	path(tiles) // Non overlapping genome bins
	tuple val(sample), val(sample_with_well_coordinate), path(allCells), path(metCG)
output:
	tuple val(sample_with_well_coordinate), path("${sample_with_well_coordinate}.CG.score.mtx"), path("${sample_with_well_coordinate}.barcodes.tsv"), path("${sample_with_well_coordinate}.CG.features.tsv"), emit: mtx, optional: true
tag "$sample_with_well_coordinate"
label 'process_single'
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
	tuple val(sample_with_well_coordinate), path("${sample_with_well_coordinate}.CH.mtx"), path("${sample_with_well_coordinate}.barcodes.tsv"), path("${sample_with_well_coordinate}.CH.features.tsv"), emit: mtx, optional: true
tag "$sample_with_well_coordinate"
label 'process_medium'
script:
"""
	create_mtx.py --bedfile $tiles --met_calls $metCH --all_cells $allCells --sample $sample_with_well_coordinate --CH
"""
}

// Generate merged CG matrix with all barcodes
process MergeMtxCG {
input:
	tuple val(sample), path(mtx), path(barcodes), path(features), path(allCells)
output:
	tuple val(sample), path("${sample}.CG.score.mtx.gz")
	tuple val(sample), path("${sample}.barcodes.tsv")
	tuple val(sample), path("${sample}.CG.features.tsv")
tag "$sample"
publishDir file(params.outDir) / "samples" / "genome_bin_matrix", pattern: "*.{mtx.gz,tsv}", mode: 'copy'
label 'process_medium_memory'
script:
"""
	merge_mtx.py --sample ${sample} --barcodes ${barcodes} --features ${features[0]} --mtx $mtx --all_cells $allCells
"""
}

// Generate merged CH matrix with all barcodes
process MergeMtxCH {
input:
	tuple val(sample), path(mtx), path(barcodes), path(features), path(allCells)
output:
	tuple val(sample), path("${sample}.CH.mtx.gz")
	tuple val(sample), path("${sample}.barcodes.tsv")
	tuple val(sample), path("${sample}.CH.features.tsv")
tag "$sample"
publishDir file(params.outDir) / "samples" / "genome_bin_matrix", pattern: "*.{mtx.gz,tsv}", mode: 'copy'
label 'process_medium_memory'
script:
"""
	merge_mtx.py --sample ${sample} --barcodes ${barcodes} --features ${features[0]} --mtx $mtx --all_cells $allCells
"""
}

// Create .allc files for allcools from internal methylation coverage parquet files and tabix index
process CreateALLC {
input:
	tuple val(sample), val(sample_with_well_coordinate), path(allCells), path(met)
output:
	tuple val(sample), path("{CG,CH}/*.allc.tsv.gz"), optional: true
	tuple val(sample), path("{CG,CH}/*.tbi"), optional: true
tag "$sample_with_well_coordinate"
publishDir { outDir }, pattern: "{CG,CH}/*.allc.tsv.gz", mode: 'copy'
publishDir { outDir }, pattern: "{CG,CH}/*.tbi", mode: 'copy'
label 'process_medium_memory'
script:
    outDir = file(params.outDir) / "samples" / "methylation_coverage" / "allc" / sample.tokenize('.')[0]
"""
	write_allc.py --met_calls $met --all_cells $allCells --sample $sample_with_well_coordinate
	for file in {CG,CH}/*.tsv.gz; do 
		if [[ -f "\${file}" ]]; then 
		tabix -@ ${task.cpus} -b2 -e2 -s1 \${file};
		fi
	done
"""
}

process CreateBismark {
input:
	tuple val(sample), val(sample_with_well_coordinate), path(allCells), path(met)
output:
	tuple val(sample), path("{CG,CH}/*.cov.gz"), optional: true
tag "$sample_with_well_coordinate"
publishDir { outDir }, pattern: "{CG,CH}/*.cov.gz", mode: 'copy'
label 'process_medium_memory'
script:
    outDir = file(params.outDir) / "samples" / "methylation_coverage" / "cov" / sample.tokenize('.')[0]
"""
	write_bismark.py --met_calls $met --all_cells $allCells --sample $sample_with_well_coordinate
"""
}

process CreateAmethyst {
input:
	tuple val(sample), val(sample_with_well_coordinate), path(allCells), path(metCG), path(metCH)
output:
	tuple val(sample), path("*_cov.h5"), optional: true // test datasets may have no passing cells for one well
tag "$sample_with_well_coordinate"
publishDir { outDir }, pattern: "*_cov.h5", mode: 'copy'
label 'process_medium_memory'
script:
    outDir = file(params.outDir) / "samples" / "methylation_coverage" / "amethyst" / sample.tokenize('.')[0]
"""
	write_amethyst.py --met_cg $metCG --met_ch $metCH --all_cells $allCells --sample $sample_with_well_coordinate
"""
}


workflow MATRIX_GENERATION {
take:
    allCells // Information about all cell barcodes
    metCG // Per TN5 CG methylation stats
    metCH // Per TN5 CH methylation stats
    genome // Reference genome

main:
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

    if (params.matrixGenerationCG) {
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
            .join(allCells) // Output matrix with barcodes in same order as metadata file
            .dump(tag: 'mtxCGPerSample')
            .set { mtxCGPerSample }
        MergeMtxCG(mtxCGPerSample)
    }
    
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

    if (params.matrixGenerationCH) {
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
            .join(allCells) // Output matrix with barcodes in same order as metadata file
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
    
    if (params.amethystOut) {
        // join ch and cg channels for CreateAmethyst
        cgPass
            .join(chPass, by: [0, 1, 2], failOnDuplicate: true, failOnMismatch: true)
            .dump(tag: 'metPass')
            .set { metPass }
        CreateAmethyst(metPass)
    }
    
}
