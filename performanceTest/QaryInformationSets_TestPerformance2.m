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
procedure TestPerformanceC_x(fileOutput, C, typeAlgMethod)

    kernel := C`Kernel;
    cosetRep := C`CosetRepresentatives;

    maxDimPartialK := Dimension(kernel);
    minDimPartialK := 0;

    PrintFile(fileOutput, "Partial kernel dimensions tested:");
    PrintFile(fileOutput, [minDimPartialK..maxDimPartialK]);

    timePerformance := [];

    for dimPartialK := maxDimPartialK to minDimPartialK by -1 do

        //partial kernel and coset representatives
        partialK := Subcode(kernel, dimPartialK);
        compK := CodeComplement(kernel, partialK);
        partialRep := [k + c : k in Set(compK), c in cosetRep ];
        newCode := QaryCode(partialK, partialRep : IsFinalKernel := true);
        case typeAlgMethod:
            when "BruteForce1":
                tstart := Cputime();
                for iteration in [1..20] do
                    OutputAllInformationSets := AllInformationSetsV2(newCode);
                end for;
                tend := Cputime(tstart);
            when "BruteForce2":
                tstart := Cputime();
                for iteration in [1..20] do
                    OutputAllInformationSets := AllInformationSetsV3(newCode);
                end for;
                tend := Cputime(tstart);
            when "CosetsRep":
                tstart := Cputime();
                for iteration in [1..20] do
                    OutputAllInformationSets := AllInformationSetsV1(newCode);
                end for;
                tend := Cputime(tstart);
        end case;

        Append(~timePerformance, tend/20);
        output := Sprintf("Partial kernel dimension = %o, time (seconds) = %o ", 
                        dimPartialK, tend/20);
        PrintFile(fileOutput, output);

    end for;

    PrintFile(fileOutput, timePerformance);
    
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


print "test 11: q=2, n=16, kerDim=4, #C=2048, #rep=128, nonlinear, 0 notin C";
load "./data/CHamq2n16ker4t_seqgen.m";
L := CHamq2n16ker4t_seq;
C := QaryCode(L);

fileOutput1 := createFile("TestPerformanceAllInformationSets");
PrintFile(fileOutput1, C);

PrintFile(fileOutput1, "Test IsInfSet method Coset Representatives V1");
TestPerformanceC_x(fileOutput1, C,  "CosetsRep");

PrintFile(fileOutput1, "Test IsInfSet method BruteForce1 V2");
TestPerformanceC_x(fileOutput1, C,  "BruteForce1");

PrintFile(fileOutput1, "Test IsInfSet method BruteForce2 V3");
TestPerformanceC_x(fileOutput1, C,  "BruteForce2"); 




