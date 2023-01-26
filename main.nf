#!/usr/bin/env nextflow

// Declare syntax version
nextflow.enable.dsl=2

include { HCM } from './workflows/HCM.nf'

workflow {
    HCM()
}