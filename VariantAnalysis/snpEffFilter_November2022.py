import sys

# Set input and output variables
input = sys.argv[1]
orthologyFile = open(sys.argv[2],"r")
goTermFile = open(sys.argv[3],"r")
inHouseFile = open(sys.argv[4],"r")
# Set various options: keep low, keep modifier. Default
# is off, putting the options in the command line will
# automatically turn to true, no further specification
# needed.
keepLow = False
keepMod = False
index = 0
while index < len(sys.argv):
    if sys.argv[index] == "-low":
        keepLow = True
    if sys.argv[index] == "-mod":
        keepMod = True
    index += 1
# Begin iterating through input file. First, must collect
# annotation information
with open(input,"r") as inpFile:
    for line in inpFile:
        removeFlag = False
        if "#" in line:
            continue
        # Split line by tab value. Annotation is in index 7
        lineList = line.split('\t')
        lineList_outputSegment = lineList[0]+'\t'+lineList[1]+'\t'+lineList[3]+'\t'+lineList[4]
        for houseLine in inHouseFile:
            houseList = houseLine.split('\t')
            houseString = houseList[0] + '\t' + houseList[1] + '\t' + houseList[2] + '\t' + houseList[3]
            if houseString == lineList_outputSegment:
                removeFlag = True
                break
        inHouseFile.seek(0)
        if removeFlag:
            continue
        annotations = lineList[7]
        # Split annotation by "|"
        annotationList = annotations.split("|")
        # Filter annotationList based on the options provided
        # This code will increment through the annotations, find
        # the VEP (HIGH or MODERATE) and collect relevant information
        # based on its relative position to the VEP
        annotationIndex = 0
        filteredAnnotations = []
        while annotationIndex < len(annotationList):
            if annotationList[annotationIndex] == "HIGH" or annotationList[annotationIndex] == "MODERATE":
                buildAnnotation = []
                effect = annotationList[annotationIndex - 1]
                VEP = annotationList[annotationIndex]
                wbID = annotationList[annotationIndex + 2]
                change = annotationList[annotationIndex + 7]
                buildAnnotation.append(effect)
                buildAnnotation.append(VEP)
                buildAnnotation.append(wbID)
                buildAnnotation.append(change)
                filteredAnnotations.append(buildAnnotation)
            if keepLow:
                if annotationList[annotationIndex] == "LOW":
                    buildAnnotation = []
                    effect = annotationList[annotationIndex - 1]
                    VEP = annotationList[annotationIndex]
                    wbID = annotationList[annotationIndex + 2]
                    change = annotationList[annotationIndex + 7]
                    buildAnnotation.append(effect)
                    buildAnnotation.append(VEP)
                    buildAnnotation.append(wbID)
                    buildAnnotation.append(change)
                    filteredAnnotations.append(buildAnnotation)
            if keepMod:
                if annotationList[annotationIndex] == "MODIFIER":
                    buildAnnotation = []
                    effect = annotationList[annotationIndex - 1]
                    VEP = annotationList[annotationIndex]
                    wbID = annotationList[annotationIndex + 2]
                    change = annotationList[annotationIndex + 7]
                    buildAnnotation.append(effect)
                    buildAnnotation.append(VEP)
                    buildAnnotation.append(wbID)
                    buildAnnotation.append(change)
                    filteredAnnotations.append(buildAnnotation)
            annotationIndex += 1
        # Collect annotation info: effect,VEP,wbID,change
        # effect is index 0, VEP is index 1, wbID is index
        # 2, change is index 3
        annotationIndex = 0
        effect = ""
        VEP = ""
        wbID = ""
        change = ""
        geneName = ""
        locusName = ""
        humanGenes = ""
        omim = ""
        goDesc = ""
        while annotationIndex < len(filteredAnnotations):
            annotation_outputSegment = ""
            annotation_info = filteredAnnotations[annotationIndex]
            effect = annotation_info[0]
            VEP = annotation_info[1]
            wbID = annotation_info[2]
            change = annotation_info[3]
            # Search up info in orthology and go term files
            for orthoLine in orthologyFile:
                orthoLineList = orthoLine.split('\t')
                if wbID == orthoLineList[0]:
                    geneName = orthoLineList[1]
                    locusName = orthoLineList[2]
                    humanGenes = orthoLineList[4]
                    omim = orthoLineList[5]
                    omim = omim.rstrip()
                    break
            orthologyFile.seek(0)
            for termLine in goTermFile:
                termLineList = termLine.split('\t')
                if wbID == termLineList[0]:
                    goDesc = termLineList[1] + '\t' + termLineList[2] + '\t' + termLineList[3]
                    goDesc = goDesc.rstrip()
                    break
            goTermFile.seek(0)
            # Format output line        
            annotation_outputSegment += '\t' + effect + '\t' + VEP + '\t' + wbID + '\t' + change + '\t' + locusName + '\t' + geneName + '\t' + humanGenes + '\t' + goDesc + '\t' + omim
            # Format output line
            outputLine = lineList_outputSegment+annotation_outputSegment+'\n'
            print(outputLine)
            annotationIndex += 1
orthologyFile.close()
goTermFile.close()
