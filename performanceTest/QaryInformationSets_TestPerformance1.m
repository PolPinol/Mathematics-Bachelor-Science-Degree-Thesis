/****************************************************************/
/*                                                              */
/* Project name: Qary nonlinear codes in MAGMA	                */
/* Test file name: QaryInformationSets_TestPerformance.m        */
/*                                                          	*/
/* Comments: Performance test for the function ...              */
/*           included in the QaryCodes_InformationSets.m file   */
/*                                                              */
/* Authors: P. Piñol and M. Villanueva                          */
/*                                                          	*/
/* Revision version and last date: 1.0,  2019/06/20         	*/
/*                                 1.1,  2019/09/19             */
/*                                                              */
/****************************************************************/

/****************************************************************/
/* Performance test description:                                */
/*                                                              */
/*   The output file returns the time in seconds of running     */
/*   function C^x using the extended cosets.                    */                            
/*   It also return the time in seconds to obtain the same      */
/*   result with the same function but using brute force.       */
/*                                                              */ 
/*   These results determine up to which dimension of the kernel*/
/*   it is more efficient to use extended cosets instead of     */
/*   a brute force method.                                      */
/*                                                              */ 
/****************************************************************/

Attach("QaryCodes_Core.m");
Attach("QaryCodes_Extension.m");
Attach("QaryCodes_Distances.m");
Attach("QaryCodes_InformationSets.m");

/****************************************************************/
/*                                                              */
/* Function name: createFile                                    */
/* Parameters: fileOutput                                       */
/* Function description: Create a magma file named fileOutput,  */
/*   where system server information is printed.                */
/* Input parameters description:                                */
/*   - fileOutput: string                                       */
/* Output parameters description:                               */
/*   - the magma file                                           */
/*                                                              */
/****************************************************************/
function createFile(fileOutput)

    //Host information
    f := POpen("cat /proc/cpuinfo | grep model", "r");
    model_name := Gets(f);
    PrintFile(fileOutput, model_name);
    
    //MAGMA information
    memory := GetMemoryUsage();
    V, a, b := GetVersion();
    PrintFile(fileOutput, Sprintf("Magma Version: %o.%o.%o", V, a, b));
    PrintFile(fileOutput, Sprintf("Memory usage (bytes): %o", memory));

    return fileOutput;

end function;

/****************************************************************/
/*                                                              */
/* Function name: createCode                                    */
/* Parameters: recordC                                          */
/* Function description: Create a q-ary code from a record      */
/*   with a base of the kernel and the coset representatives.   */
/* Input parameters description:                                */
/*   - recordC: A record with information of the code           */
/* Output parameters description:                               */
/*   - A q-ary code                                             */
/*                                                              */
/****************************************************************/
function createCode(recordC)

    kernel := LinearCode(Matrix(recordC`Kernel));
    cosetsRepresentatives := recordC`CosetRepresentatives;

    return QaryCode(kernel, cosetsRepresentatives : IsFinalKernel := true);

end function;

/****************************************************************/
/*                                                              */
/* Function name: TestPerformanceC_x                            */
/* Parameters: fileOutput, C, permutations, typeAlgMethod       */
/* Procedure description: Construct the test performance        */
/*   of the function C^x considering artifical codes with       */
/*   diferent dimension of the kernel from the same code.       */                                            
/* Input parameters description:                                */
/*   - fileOutput: A file                                       */
/*   - C: A q-ary code                                          */
/*   - permutations: A sequence of permutations                 */
/*   - typeAlgMethod: Method to test the performance            */
/*                                                              */
/****************************************************************/
procedure TestPerformanceC_x(fileOutput, C, informationSets, typeAlgMethod)

    OutputList := [IsInformationSet(C,I) : I in informationSets];
    kernel := C`Kernel;
    cosetRep := C`CosetRepresentatives;

    maxDimPartialK := Dimension(kernel);
    minDimPartialK := 0;

    PrintFile(fileOutput, "Partial kernel dimensions tested:");
    PrintFile(fileOutput, [minDimPartialK..maxDimPartialK]);

    i := 1;
    for I in informationSets do
    
        PrintFile(fileOutput, "Results for the information set");
        PrintFile(fileOutput, I);
        timePerformance := [];

        for dimPartialK := maxDimPartialK to minDimPartialK by -1 do

            //partial kernel and coset representatives
            partialK := Subcode(kernel, dimPartialK);
            compK := CodeComplement(kernel, partialK);
            partialRep := [k + c : k in Set(compK), c in cosetRep ];
            newCode := QaryCode(partialK, partialRep : IsFinalKernel := true);
            case typeAlgMethod:
                when "BruteForce3":
                    tstart := Cputime();
                    for iteration in [1..20] do
                        OutputIsInformationSet := IsInformationSetBF3(newCode,I);
                    end for;
                    tend := Cputime(tstart);
            end case;

            Append(~timePerformance, tend/20);
            output := Sprintf("Partial kernel dimension = %o, time (seconds) = %o ", 
                            dimPartialK, tend/20);
            PrintFile(fileOutput, output);

            assert OutputIsInformationSet eq OutputList[i];

        end for;

        PrintFile(fileOutput, timePerformance);
        i := i+1;

    end for;

end procedure;

/****************************************************************/
/* Code Cq3n13ker8a: qaryrec record                             */
/*   - Length: 13                                               */
/*   - Minimum Distance: 3                                      */
/*   - Number of Codewords:  59049                              */
/*   - IsLinear: false                                          */
/*   - Kernel dimension: 8                                      */
/*   - #Coset representatives: 9                                */
/*                                                              */
/* Comments:  Hamming code over GF(3)                           */
/*            Test #12 from QaryGroupAction_BB_test.m           */
/*                                                              */
/****************************************************************/

print "test 1: Universe code of length 10 over GF(3)";
universeCode := UniverseCode(GF(3), 10);
C := QaryCode(universeCode);
fileOutput1 := createFile("TestPerformanceIsInformationSet");
PrintFile(fileOutput1, C);

I1 := [1..Length(C)];
infoSets := [I1];
PrintFile(fileOutput1, "Test IsInfSet method BruteForce3");
TestPerformanceC_x(fileOutput1, C, infoSets, "BruteForce3");


print "test 2: Repetition code of length 12 over GF(4)";
repetitionCode := RepetitionCode(GF(4), 12);
C := QaryCode(repetitionCode);
fileOutput1 := createFile("TestPerformanceIsInformationSet");
PrintFile(fileOutput1, C);

I1 := [1];
infoSets := [I1];
PrintFile(fileOutput1, "Test IsInfSet method BruteForce3");
TestPerformanceC_x(fileOutput1, C, infoSets, "BruteForce3");


print "test 3: Zero code of length 13 over GF(8)";
zeroCode := ZeroCode(GF(8), 13);
C := QaryCode(zeroCode);
fileOutput1 := createFile("TestPerformanceIsInformationSet");
PrintFile(fileOutput1, C);
I1 := [];
infoSets := [I1];
PrintFile(fileOutput1, "Test IsInfSet method BruteForce3");
TestPerformanceC_x(fileOutput1, C, infoSets, "BruteForce3");


print "test 4: q=3, n=13, kerDim=8, #C=59049, #rep=9, nonlinear"; 
load "./data/Cq3n13ker8c_seqgen.m";
L := Cq3n13ker8c_seq;
C := QaryCode(L);
fileOutput1 := createFile("TestPerformanceIsInformationSet");
PrintFile(fileOutput1, C);
I1 := [1,2,3,4,5,6,7,8,10,11];
infoSets := [I1];
PrintFile(fileOutput1, "Test IsInfSet method BruteForce3");
TestPerformanceC_x(fileOutput1, C, infoSets, "BruteForce3");

print "test 5: q=3, n=13, kerDim=9, #C=59049, #rep=3, nonlinear"; 
load "./data/Cq3n13ker9a_seqgen.m";
L := Cq3n13ker9a_seq;
C := QaryCode(L);
fileOutput1 := createFile("TestPerformanceIsInformationSet");
PrintFile(fileOutput1, C);
I1 := [1,2,3,4,5,6,7,8,10,11];
infoSets := [I1];
PrintFile(fileOutput1, "Test IsInfSet method BruteForce3");
TestPerformanceC_x(fileOutput1, C, infoSets, "BruteForce3");

print "test 6: q=4, n=16, kerDim=1, #C=64, #rep=16, nonlinear"; 
load "./data/CHadq4n16ker1b_seqgen.m";
L := CHadq4n16ker1b_seq;
C := QaryCode(L);
fileOutput1 := createFile("TestPerformanceIsInformationSet");
PrintFile(fileOutput1, C);
I1 := [2,7,8];
infoSets := [I1];
PrintFile(fileOutput1, "Test IsInfSet method BruteForce3");
TestPerformanceC_x(fileOutput1, C, infoSets, "BruteForce3");

