#!/usr/bin/env nextflow

// Declare syntax version
nextflow.enable.dsl=2

include { HCM_cardiac_enhs } from './workflows/HCM_cardiac_enhs.nf'

workflow {
    HCM_cardiac_enhs()
}