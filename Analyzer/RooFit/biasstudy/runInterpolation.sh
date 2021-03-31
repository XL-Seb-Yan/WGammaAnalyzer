#!/bin/bash
root -l -q makeRooMultiPdfWorkspaceSigW.C
root -l -q plotSignalPDF.C
rm signal_pdfs_*
