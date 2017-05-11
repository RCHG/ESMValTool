'''
This script is part of the ESMValTool distribution.  It's been added as part of the incorporation
of the Sphinx documentation generator.  Sphinx was originally developed for documenting Python code, 
and one of its features is that it is able - using the so-called autodoc extension - to extract 
documentation strings from Python source files and use them in the documentation it generates.  

The autodoc feature apparently does not exist for NCL source files (such as those which are used in 
ESMValTool), but it has been mimicked (or - more-or-less - reverse-engineered) here via this script,
which walks through a subset of the ESMValTool NCL scripts, extracts function names, argument lists 
and descriptions (from the comments immediately following the function definition), and assembles 
them in a subdirectory of doc/sphinx/source.  These output files are in the so-called reStructuredText 
format (see, e.g., http://docutils.sourceforge.net/rst.html), which is the markup language used by 
Sphinx; running make in doc/sphinx builds the ESMValTool documentation from them, as noted above. 

Created on July 14, 2015

@author: jeremy.walton@metoffice.gov.uk
'''

import os
import glob
import re
import string
import collections
import sys
import datetime


def openFile(fileName, mode):
    
    ''' Opens a file fileName using mode (reading or writing) and return the pointer if successful.  '''
    
    # Try to open the file.
    try:
        filePointer = open(fileName, mode)
        return filePointer
    except IOError:
        print "Couldn't open", fileName
        return None


def writeTimeStamp(filePointer):

    ''' Write a comment to the file (pointed to by filePointer) to show it's been
    automatically generated at this time.  '''
    
    now = datetime.datetime.now()
    filePointer.write('.. This file has been automatically generated by ' + __file__ + \
          ' on ' + now.strftime("%Y-%m-%d %H:%M:%S") +'\n\n')



def makeParamDetails(params):
    
    ''' Create a list of parameter names and types from the params string.  '''
    
    # We'll store the parameter names and their types in a dictionary.  Note that it has to be 
    # an ordered dictionary, because later on we want to pull the entries out in the same order 
    # that we added them.   
    paramDetails = collections.OrderedDict()
    for param in params:
        
        # Extract the type if it's specified, otherwise default to integer (say).
        if ':' in param:
            [pname, ptype] = param.split(':')
        else:
            pname = param
            ptype = 'integer'
        
        # If the parameter is an array, we only want its name in the description.
        pname = pname.split('[')[0]
        pname = pname.strip()
        
        # Tie the name and the type of the parameter together.
        paramDetails[pname] = ptype
        
    return paramDetails


def processParams(params, inp, oup):
    
    ''' Extract the parameter names and types from the params string, pull their 
    descriptions out from the input file and reformat the lot in the output.  '''

    # Get the names and types.
    paramDetails = makeParamDetails(params)

    # We assume we're at the line before the first parameter description.  Bump it, then check to see
    # if we're really at the right location and issue a warning if not.
    line = inp.next()
    if paramDetails.keys()[0] not in line:
        print "Warning - parameter " + paramDetails.keys()[0] + " not found in this line:\n" + line
    
    # Want a blank line just before parameter descriptions.
    oup.write('\n')   
       
    # Loop over all parameters in the argument list.  
    for i, pname in enumerate(paramDetails.keys()):
        
        # Now assemble the description from the line(s).  We assume that there's a description for each 
        # parameter, and that the descriptions are in the same order as the parameters in the function
        # declaration.  The descriptions could be spread over more than one line.  We have some 
        # very rudimentary checking below to detect deviations from these assumptions (or requirements). 
        if pname in line:
            
            # Get the text in the line which follows the first occurrence (reading from the left) 
            # of the parameter name, then strip trailing spaces (including the CR).
            pdesc = line.split(pname, 1)[1]
            pdesc = pdesc.rstrip()
                                          
            # The description could continue on the following lines, which need to be concatenated
            # together.  For all except the last parameter, the end of the description is signaled 
            # by the name of the next parameter.  For the last (or maybe the only) parameter, it's 
            # signaled by a blank line.  
            line = inp.next()
            if i < len(paramDetails.keys())-1:
                pnext = paramDetails.keys()[i+1]     
                if pnext not in line:
                    
                    # If the line contains a colon, signal it as unexpected, because it could be 
                    # the description of a parameter which is different from the one we're looking 
                    # for.
                    if ':' in line:
                        print "Warning - unexpected line:"
                        print line.rstrip()
                        print "encountered when searching for description for parameter " + pnext
                        
                    # Do the concatentation, stripping whitespace (including the CR) as we go.
                    while not pnext in line:
                        pdesc += " " + line.replace(';;', '  ', 1).strip()
                        line = inp.next()
            else:
                while not line.replace(';;', '  ', 1).isspace():
                    pdesc += " " + line.replace(';;', '  ', 1).strip()                 
                    line = inp.next()

            # Ensure the description starts with a colon.
            if pdesc[0] != ':':
                pdesc = ':' + pdesc             
            
            # Write out the complete description of this parameter. 
            oup.write('   :param ' + paramDetails[pname] + ' ' + pname + pdesc + '\n') 
             
        else:
            print "Warning - parameter " + paramDetails.keys()[i] + " not found in this line:\n" + line
                            
    # Want a blank line just after parameter descriptions.
    oup.write('\n')   


def findArgument(inp):
    
    ''' Find the start of the Arguments list.  '''
    
    line = inp.next()
    count = 1
    while not 'Arguments' in line:
        line = inp.next()
        
        # We assume we're going to find this within two lines of the original location of the input
        # - stop looking if we don't.
        count = count + 1
        if count > 2:
            return False
    
    return True


def parseFile(inFileName, outFileName):
    
    ''' Processes an ncl file and produces an rst file as output, which contains documentation of the 
    ncl functions in a form suitable for input to the Sphinx documentation generator.  ''' 
    
    # Open the files.
    inp = openFile(inFileName, "r")
    if inp is None:
        return

    oup = openFile(outFileName, "w")
    if oup is None:
        return
    
    # Write a comment to the file to show it's been automatically generated.
    writeTimeStamp(oup)

    # We assume the file name has the form /path/to/foo.ncl, and the module name is foo.  Pull it out,
    # and write it to the output file as the title.
    modName = os.path.splitext(os.path.basename(inFileName))[0]
    
    oup.write(':mod:' + '`' + modName + '`' + '\n')
    oup.write("=" * (7+len(modName)) + '\n')
    
    for line in inp:
        
        # Is this the start of a function?
        if re.match('^function', line) or re.match('^procedure', line):
            
            # The function could have parameters on the following lines.  Concatenate them
            # up until the closing bracket, stripping whitespace (including the CR) as we go.
            fname = line.rstrip()
            while ')' not in fname:
                line = inp.next()
                fname += " " + line.strip()

            # Some ncl files have backslashes in the function declaration to indicate continuation to the 
            # next line (even though this isn't necessary in ncl).  These will mess up our processing of 
            # the argument list, and don't look good in the doc. so we pull them out here.
            fname = string.replace(fname, '\\', '')

            # Write the line out from the word 'function' onwards, and suitably decorated for rst.  
            # Need the CR at the end, as we've been pulling that off throughout the assembly of this line.  
            oup.write('.. function:: ' + fname[len('function')+1:] + '\n')
            
            # Print a message to say we're processing this function.
            print "Extracting details from " + fname.split('(')[0]
            
            # Now extract the list of parameters from the function declaration.  First, pull 
            # out the text between the brackets, then split that into individual parameter names.
            plist = fname.split('(')[1].split(')')[0]
            params = plist.split(',')
                               
            # Position the input just after the line containing 'Arguments'.  
            if not findArgument(inp):
                print "Warning - argument list not found for " + fname
            else:
                
                # Here's where we check whether this function has any parameters.  If it doesn't,
                # then we don't need to process any.               
                if len(plist) > 0:
                    # Read the parameter descriptions and reformat them before writing them out.
                    line = processParams(params, inp, oup)
                
                # We assume the first batch of comments immediately following the function are 
                # part of the documentation.
                line = inp.next()
                while re.match('^;;', line):
                    
                    # Make the line to be output by replacing comments with spaces.
                    outLine = string.replace(line, ';;', '  ', 1)
                    
                    # If this is the start of the modification history, then add an extra
                    # colon and line feed, which signals that the following text - i.e. the 
                    # list of modifications - is preformatted (so that, e.g. line breaks are
                    # preserved).
                    if 'Modification history' in outLine:
                        # First, pull off the trailing whitespace (including CR).
                        outLine = outLine.rstrip()
                        
                        # Some of the headings already end in a colon, some don't.  
                        if ':' in outLine:
                            outLine = outLine + ':'
                        else:
                            outLine = outLine + '::'
                        
                        # Now replace the CR, and add an extra one. 
                        outLine = outLine + '\n\n' 
                    
                    # Write out this line and get the next one.
                    oup.write(outLine)
                    line = inp.next()
                
    # Close the files.
    inp.close()
    oup.close()    

def processFiles():
    
    ''' Does all the processing for the files.  (This is all in a single routine for ease of switching 
    full functionality on and off when testing.) '''
    
    # Do some rudimentary checking of where this script is being run from, because we're going to be 
    # using relative paths below to find the directories containing the input & output.
    if not os.path.exists("../../doc/sphinx/scripts"):
        print "Error - this script should be run from within the doc/sphinx directory using the command"
        print "% python scripts/process_ncl_docs.py"
        sys.exit()

    # List the directories containing input files, then loop over them.
    inDirs = ["../../diag_scripts/lib/ncl/", "../../plot_scripts/ncl/"]
    for inDir in inDirs:

        # Get the useful part of the name of the input directory (we assume full name
        # is ../../../foo/bar, where foo is the useful part of the name).
        dirName = inDir.split('/')[2]

        # Make an index file for this directory, and open it.
        indexFile = "source/" + dirName + "/index.rst"
        oup = openFile(indexFile, "w")
        if oup is None:
            return

        # Write a comment to the file to show it's been automatically generated.
        writeTimeStamp(oup)

        # Write the preamble to the file.
        oup.write('NCL ' + dirName + ' library' + '\n')
        oup.write("=" * (12+len(dirName)) + '\n\n')
        oup.write('Contents: \n\n')
        oup.write('.. toctree:: \n')
        oup.write('   :maxdepth: 2 \n\n')

        # Form the output directory name from the input directory name.
        outDir = "source/" + dirName + '/'

        # Find all the ncl files in the input directory, and loop over them.
        inFiles = glob.glob(inDir + '*.ncl')
        print inFiles
        for nclFile in inFiles:
            print "Processing " + nclFile
            fileName = os.path.basename(nclFile).replace('.ncl', '')
            oup.write('   ' + fileName + '\n')
            rstFile = outDir + os.path.basename(nclFile).replace('.ncl', '.rst')
            parseFile(nclFile, rstFile)

        # Close the index file.
        oup.close()


if __name__ == '__main__':

    processFiles()
#    nclFile = "monsoon_panels.ncl"
#    outDir = "."
#    print "Processing " + nclFile
#    rstFile = outDir + os.path.basename(nclFile).replace('.ncl', '.rst')
#    parseFile(nclFile, rstFile)
