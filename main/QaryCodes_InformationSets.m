///////////////////////////////////////////////////////////////////////////////
/////////    Copyright 2019-20 Pol Piñol under the supervision          ///////
/////////              of Jaume Pujol and Mercè Villanueva              ///////
/////////                                                               ///////
/////////    This program is distributed under the terms of GNU         ///////
/////////               General Public License                          ///////
///////////////////////////////////////////////////////////////////////////////

/******************************************************************************
    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
******************************************************************************/

/************************************************************/
/*                                                          */
/* Project name: Qary nonlinear codes in MAGMA              */
/* File name: QaryCodes_InformationSets.m                   */
/*                                                          */
/* Comment: Package developed within the CCSG group         */
/*                                                          */
/* Revision version and last date: version 1.0 08/06/2020   */
/*                                                          */
/************************************************************/
//Uncomment freeze when package finished
//freeze;

/* PACKAGE VERSION */
intrinsic QaryCodes_InformationSets_version() -> SeqEnum
{Return the current version of this package.}
    version := [1,0];
    return version;
end intrinsic;

//import "QaryCodes_Core.m": NewCodeFld;

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////               INFORMATION SETS FUNCTIONS                        ///////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

ParametersInfoSet := function(sizeCode, finiteField)
    p := Characteristic(finiteField);
    s := Degree(finiteField);
    t, e := Valuation(sizeCode, p);
    t := Floor(t/s);

    return t, e;
end function;

DiagonalizeKernelFromInfoSetK := function(K, IK)   
    M := GeneratorMatrix(K);
    k := Dimension(K);
    for j in [1..#IK] do

		columnPivot := IK[j];
		rowPivot := j;	
        	            
		// search for nonzero element in columnPivot column
        while (rowPivot lt k) and (M[rowPivot][columnPivot] eq 0) do
     	    rowPivot := rowPivot+1;
        end while;
	    // swap rows if nonzero element found
	    SwapRows(~M,j,rowPivot);
            
        // diagonalize columnPivot column
        M[j] := M[j]/M[j][IK[j]];
        for i := j-1 to 1 by -1 do
            if (M[i][IK[j]] ne 0) then
   	            M[i] := M[i] - M[i][IK[j]] * M[j];
            end if;
        end for;
        for i := j+1 to k do
            if (M[i][IK[j]] ne 0) then
	            M[i] := M[i] - M[i][IK[j]] * M[j];
            end if;
        end for;
    end for;
    
    return M;
end function;


/****************************************************************/
/*                                                              */
/* Function name: IsSystematic                                  */
/* Parameters: C                                                */
/* Function description:  Return true if and only if the        */
/* q-ary code C of length n is systematic, that is, there is    */
/* a set of $t$ coordinate positions I C {1,..,n}, called       */
/* information set, such that |C|=|C_I|=q^t, where              */
/* C_I={ v_I : v in C } and v_I denote the restriction of       */
/* the vector v to the coordinates in I.                        */
/* Input parameters description:                                */
/*   - C: A q-ary code                                          */
/* Output parameters description:                               */
/*   - true if I is Systematic, false otherwise                 */
/*   - An information set I for C                               */
/*   - A tuple < I_k, I_r > where I_k is an information         */
/*     set for the kernel K(C) and I=I_k U I_r is an            */
/*     information set for C.                                   */
/*                                                              */
/* Signature: (<CodeFld> C) -> BoolElt, SeqEnum, Tup            */
/*                                                              */
/****************************************************************/
intrinsic IsSystematic(C::CodeFld) -> BoolElt, SeqEnum, Tup
{Return true if and only if the q-ary code C of length n is systematic, that is, there is a set of t coordinate positions I C [1,..,n], called information set, such that |C|=|C_I|=q^t, where C_I=[ v_I : v in C ] and v_I denote the restriction of the vector v to the coordinates in I.}

    if(C`IsLinear) then
        IK := InformationSet(C`Kernel);
        return true, IK, <IK,[]>;
    end if;
    
    n := Length(C);
    t, e := ParametersInfoSet(#C, C`BaseField);
    k :=  Dimension(C`Kernel);
    r := t-k;
    V := InformationSpace(C`Kernel);
    numCosetRep := #C`CosetRepresentatives;
    
    // #C is not a power of q 
    if e ne 1 then
        return false, [], <[],[]>;
    end if;

    // #C is 1, so C is the zero code
    if t eq 0 then 
        return true, [], <[], []>;
    end if;

    if k ne 0 then
        for IK in AllInformationSets(C`Kernel) do
            GIK := DiagonalizeKernelFromInfoSetK(C`Kernel, IK);
            newCosetRepresentatives := [v - V![v[i] : i in IK] * GIK : v in C`CosetRepresentatives];
            for IR in Subsets({1..n} diff Set(IK), r) do
                setRepresentatives := { [v[i] : i in IR] : v in newCosetRepresentatives };
                if #setRepresentatives eq numCosetRep then
                    return true, Sort(IK cat Setseq(IR)), <IK, Setseq(IR)>;
                end if; 
            end for; 
        end for;
    else
        for I in Subsets({1..n}, t) do
            CIcodewords := { [v[i] : i in I] : v in C`CosetRepresentatives };
            if #CIcodewords eq #C then
                newI := Setseq(I);
                return true, newI, <[], newI>;
            end if;
        end for;
    end if;
      
    return false, [], <[],[]>;

end intrinsic;




/****************************************************************/
/*                                                              */
/* Function name: IsInformationSet                              */
/* Parameters: C, I                                             */
/* Function description: Return true if and only if the set     */
/* I C {1,...,n} of coordinate positions is an information set  */
/* for the q-ary code C of length n. An information set         */
/* for C is an ordered set of t coordinate positions            */
/* I C {1,...,n} such that |C|=|C_I|=q^t, where                 */
/* C_I={ v_I | v in C } and v_I denote the restriction of the   */
/* vector v to the coordinates in I.                            */
/* The set I can be given as a set or sequence of integers,     */
/* or as a tuple < I_k, I_r > with two sets or sequences        */
/* of integers. In the case that I is given as a tuple, the     */
/* function returns true if and only if I_k is an information   */
/* set for the kernel K(C) and I=I_k U I_r is an information    */
/* set for C.                                                   */
/*   - C: A q-ary code                                          */
/*   - I: A subset of coordinate possitions                     */
/* Output parameters description:                               */
/*   - true if I is an information set for C, false otherwise   */
/*                                                              */
/* Signature: (<CodeFld> C, <SeqEnum> I) -> BoolElt             */
/*                                                              */
/****************************************************************/
intrinsic IsInformationSet(C::CodeFld, I::SeqEnum) -> BoolElt
{Return true if and only if the set I C [1,...,n] of coordinate positions is an information set for the q-ary code C of length n. An information set for C is an ordered set of t coordinate positions I C [1,...,n] such that |C|=|C_I|=q^t, where C_I=[ v_I | v in C ] and v_I denote the restriction of the vector v to the coordinates in I. The set I can be given as a set or sequence of integers, or as a tuple < I_k, I_r > with two sets or sequences of integers. In the case that I is given as a tuple, the function returns true if and only if I_k is an information set for the kernel K(C) and I=I_k U I_r is an information set for C.}

    n := Length(C);
    require (I subset [1..n]): "Argument 2 is not a subset of coordinates";
    t, e := ParametersInfoSet(#C, C`BaseField);
    k := Dimension(C`Kernel);  
    I := Set(I);

    if(C`IsLinear) then
        if #I eq k then
            return Rank(Submatrix(GeneratorMatrix(C`Kernel), [1..k], Setseq(I))) eq k;
        end if;
        return false;
    end if;

    r := t-k;
    V := InformationSpace(C`Kernel);
    
    numCosetRep := #C`CosetRepresentatives;

    // #C is not a power of q 
    if e ne 1 then
        return false;
    // Since #C=q^t, #I has to be equal to t 
    elif t ne #I then
        return false;
    else
        if k ne 0 then
            for IK in AllInformationSets(C`Kernel) do
                if IK subset I then
                    GIK := DiagonalizeKernelFromInfoSetK(C`Kernel, IK);
                    newCosetRepresentatives := [v - V![v[i] : i in IK] * GIK : v in C`CosetRepresentatives];
                    for IR in Subsets(I diff Set(IK), r) do
                        setRepresentatives := { [v[i] : i in IR] : v in newCosetRepresentatives };
                        if #setRepresentatives eq numCosetRep then
                            return true;
                        end if; 
                    end for; 
                end if;
            end for;
        else
            CIcodewords := { [v[i] : i in I] : v in C`CosetRepresentatives };
            if #CIcodewords eq #C then
                return true;
            end if;
        end if;
    end if;
    return false;

end intrinsic;




/****************************************************************/
/*                                                              */
/* Function name: IsInformationSet                              */
/* Parameters: C, I                                             */
/* Function description: Return true if and only if the set     */
/* I C {1,...,n} of coordinate positions is an information set  */
/* for the q-ary code C of length n. An information set         */
/* for C is an ordered set of t coordinate positions            */
/* I C {1,...,n} such that |C|=|C_I|=q^t, where                 */
/* C_I={ v_I | v in C } and v_I denote the restriction of the   */
/* vector v to the coordinates in I.                            */
/* The set I can be given as a set or sequence of integers,     */
/* or as a tuple < I_k, I_r > with two sets or sequences        */
/* of integers. In the case that I is given as a tuple, the     */
/* function returns true if and only if I_k is an information   */
/* set for the kernel K(C) and I=I_k U I_r is an information    */
/* set for C.                                                   */
/*   - C: A q-ary code                                          */
/*   - I: A subset of coordinate possitions                     */
/* Output parameters description:                               */
/*   - true if I is an information set for C, false otherwise   */
/*                                                              */
/* Signature: (<CodeFld> C, <SetEnum> I) -> BoolElt             */
/*                                                              */
/****************************************************************/
intrinsic IsInformationSet(C::CodeFld, I::SetEnum) -> BoolElt
{Return true if and only if the set I C [1,...,n] of coordinate positions is an information set for the q-ary code C of length n. An information set for C is an ordered set of t coordinate positions I C [1,...,n] such that |C|=|C_I|=q^t, where C_I=[ v_I | v in C ] and v_I denote the restriction of the vector v to the coordinates in I. The set I can be given as a set or sequence of integers, or as a tuple < I_k, I_r > with two sets or sequences of integers. In the case that I is given as a tuple, the function returns true if and only if I_k is an information set for the kernel K(C) and I=I_k U I_r is an information set for C.}

    n := Length(C);
    require (I subset {1..n}): "Argument 2 is not a subset of coordinates";
    t, e := ParametersInfoSet(#C, C`BaseField);
    k := Dimension(C`Kernel);  

    if(C`IsLinear) then
        if #I eq k then
            return Rank(Submatrix(GeneratorMatrix(C`Kernel), [1..k], Setseq(I))) eq k;
        end if;
        return false;
    end if;

    r := t-k;
    V := InformationSpace(C`Kernel);
    numCosetRep := #C`CosetRepresentatives;

    // #C is not a power of q 
    if e ne 1 then
        return false;
    // Since #C=q^t, #I has to be equal to t 
    elif t ne #I then
        return false;
    else
        if k ne 0 then
            for IK in AllInformationSets(C`Kernel) do
                if IK subset I then
                    GIK := DiagonalizeKernelFromInfoSetK(C`Kernel, IK);
                    newCosetRepresentatives := [v - V![v[i] : i in IK] * GIK : v in C`CosetRepresentatives];
                    for IR in Subsets(I diff Set(IK), r) do
                        setRepresentatives := { [v[i] : i in IR] : v in newCosetRepresentatives };
                        if #setRepresentatives eq numCosetRep then
                            return true;
                        end if; 
                    end for; 
                end if;
            end for;
        else
            CIcodewords := { [v[i] : i in I] : v in C`CosetRepresentatives };
            if #CIcodewords eq #C then
                return true;
            end if;
        end if;
    end if;
    return false;

end intrinsic;



/****************************************************************/
/*                                                              */
/* Function name: IsInformationSet                              */
/* Parameters: C, I                                             */
/* Function description: Return true if and only if the set     */
/* I C {1,...,n} of coordinate positions is an information set  */
/* for the q-ary code C of length n. An information set         */
/* for C is an ordered set of t coordinate positions            */
/* I C {1,...,n} such that |C|=|C_I|=q^t, where                 */
/* C_I={ v_I | v in C } and v_I denote the restriction of the   */
/* vector v to the coordinates in I.                            */
/* Input parameters description:                                */
/*   - C: A linear code                                         */
/*   - I: A subset of coordinate possitions                     */
/* Output parameters description:                               */
/*   - true if I is an information set for C, false otherwise   */
/*                                                              */
/* Signature: (<CodeLinFld> C, <SeqEnum> I) -> BoolElt          */
/*                                                              */
/****************************************************************/
intrinsic IsInformationSet(C::CodeLinFld, I::SeqEnum) -> BoolElt
{Return true if and only if the set I C [1,...,n] of coordinate positions is an information set for the q-ary code C of length n. An information set for C is an ordered set of t coordinate positions I C [1,...,n] such that |C|=|C_I|=q^t, where C_I=[ v_I | v in C ] and v_I denote the restriction of the vector v to the coordinates in I.}

    require (I subset [1..Length(C)]): "Argument 2 is not a subset of coordinates";
    k := Dimension(C);
    I := Set(I);

    if #I eq k then
        return Rank(Submatrix(GeneratorMatrix(C), [1..k], Setseq(I))) eq k;
    end if;
    
    return false;   

end intrinsic;





/****************************************************************/
/*                                                              */
/* Function name: IsInformationSet                              */
/* Parameters: C, I                                             */
/* Function description: Return true if and only if the set     */
/* I C {1,...,n} of coordinate positions is an information set  */
/* for the q-ary code C of length n. An information set         */
/* for C is an ordered set of t coordinate positions            */
/* I C {1,...,n} such that |C|=|C_I|=q^t, where                 */
/* C_I={ v_I | v in C } and v_I denote the restriction of the   */
/* vector v to the coordinates in I.                            */
/* The set I can be given as a set or sequence of integers,     */
/* or as a tuple < I_k, I_r > with two sets or sequences        */
/* of integers. In the case that I is given as a tuple, the     */
/* function returns true if and only if I_k is an information   */
/* set for the kernel K(C) and I=I_k U I_r is an information    */
/* set for C.                                                   */
/* Input parameters description:                                */
/*   - C: A q-ary code                                          */
/*   - I: A tuple < I_k, I_r > with two sets or sequences       */
/*        of integers                                           */
/* Output parameters description:                               */
/*   - true if if I_k is an information set for the kernel K(C) */
/*     and I=I_k U I_r is an information set for C.             */
/*                                                              */
/* Signature: (<CodeFld> C, <Tuple> I) -> BoolElt               */
/*                                                              */
/****************************************************************/
intrinsic IsInformationSet(C::CodeFld, I::Tup) -> BoolElt
{Return true if and only if the set I C [1,...,n] of coordinate positions is an information set for the q-ary code C of length n. An information set for C is an ordered set of t coordinate positions I C [1,...,n] such that |C|=|C_I|=q^t, where C_I=[ v_I | v in C ] and v_I denote the restriction of the vector v to the coordinates in I. The set I can be given as a set or sequence of integers, or as a tuple < I_k, I_r > with two sets or sequences of integers. In the case that I is given as a tuple, the function returns true if and only if I_k is an information set for the kernel K(C) and I=I_k U I_r is an information set for C.}

    require #I eq 2: "Argument 2 must contain two components";
    IK := I[1];
    IR := I[2];
    require (Set(IK) meet Set(IR) eq {}): "Both components in argument 2 must be disjoint";
    n := Length(C); 
    require (IK subset [1..n]): "Argument 2 is not a subset of coordinates";
    require (IR subset [1..n]): "Argument 3 is not a subset of coordinates";
    k := Dimension(C`Kernel);  
    r := t-k;
    require (#IK eq k): "Argument 2 must contain ",k," coordinates in the first component";
    require (#IR eq r): "Argument 2 must contain ",r," coordinates in the second component";

    if(C`IsLinear) then
        return Rank(Submatrix(GeneratorMatrix(C`Kernel), [1..k], IK)) eq k;
    end if;

    t, e := ParametersInfoSet(#C, C`BaseField);
    q := #C`BaseField;
    V := InformationSpace(C`Kernel);

    // #C is not a power of q 
    if e ne 1 then
        return false;
    end if;

    // #C is 1, so C is the zero code
    if t eq 0 then 
        return true;
    end if;

    if k ne 0 then
        if Rank(Submatrix(GeneratorMatrix(C`Kernel), [1..k], IK)) eq k then
            GIK := DiagonalizeKernelFromInfoSetK(C`Kernel, IK);
            newCosetRepresentatives := [v - V![v[i] : i in IK] * GIK : v in C`CosetRepresentatives];
            setRepresentatives := { [v[i] : i in IR] : v in newCosetRepresentatives };
            if #setRepresentatives eq q^r then
                return true;
            end if;
        end if;
    else
        CIcodewords := { [v[i] : i in IR] : v in C`CosetRepresentatives };
        if #CIcodewords eq #C then
            return true;
        end if;
    end if;
    return false;

end intrinsic;







/****************************************************************/
/*                                                              */
/* Function name: InformationSpace                              */
/* Parameters: C                                                */
/* Function description: Given a systematic q-ary code C of     */
/* cardinality q^t, return the vector space U=F_q^t, which is   */
/* the space of information vectors for the code C. It is not   */
/* checked whether the code is systematic or not. The function  */
/* returns this vector space even if C is not systematic.       */
/* Input parameters description:                                */
/*   - C: A q-ary code                                          */
/* Output parameters description:                               */
/*   - vector space U=Zq^t                                      */
/*                                                              */
/* Signature: (<CodeFld> C) -> ModTupFld                        */
/*                                                              */
/****************************************************************/
intrinsic InformationSpace(C::CodeFld) -> ModTupFld
{Given a systematic q-ary code C of cardinality q^t, return the vector space U=F_q^t, which is the space of information vectors for the code C. It is not checked whether the code is systematic or not. The function returns this vector space even if C is not systematic.}
    
    t, e := ParametersInfoSet(#C, C`BaseField);
    require (e eq 1):"Argument 1 has not a power of q cardinality";
    q := #C`BaseField;

    return VectorSpace(GF(q), t);
    
end intrinsic;






/****************************************************************/
/*                                                              */
/* Function name: InformationSet                                */
/* Parameters: C                                                */
/* Function description: Given a systematic q-ary code C of     */
/* length n, return an information set for C. An information    */
/* set for C is an ordered set of t coordinate positions        */
/* I C {1,...,n} such that |C|=|C_I|=q^t, where                 */
/* C_I={ v_I | v \in C } and v_I denote the restriction of      */
/* the vector v to the coordinates in I. The information set I  */
/* is returned as a sequence of t integers.                     */
/* The function also returns a tuple < I_k, I_r > with two      */
/* sequences of integers such that I_k is an information set    */
/* for the kernel K(C) and I=I_k U I_r.                         */
/* It is not checked whether the code is systematic or not.     */
/* If the function does not succeed in finding an information   */
/* set, it returns an empty sequence and a tuple with two       */
/* empty sequences.                                             */
/* Input parameters description:                                */
/*   - C: A q-ary code                                          */
/* Output parameters description:                               */
/*   - An information set I for C                               */
/*   - A tuple < I_k, I_r > where I_k is an information         */
/*     set for the kernel K(C) and I=I_k U I_r is an            */
/*     information set for C.                                   */
/*                                                              */
/* Signature: (<CodeFld> C) -> SeqEnum, Tup                     */
/*                                                              */
/****************************************************************/
intrinsic InformationSet(C::CodeFld) -> SeqEnum, Tup
{Given a systematic q-ary code C of length n, return an information set for C. An information set for C is an ordered set of t coordinate positions I C [1,...,n] such that |C|=|C_I|=q^t, where C_I=[ v_I | v \in C ] and v_I denote the restriction of the vector v to the coordinates in I. The information set I is returned as a sequence of t integers. The function also returns a tuple < I_k, I_r > with two sequences of integers such that I_k is an information set for the kernel K(C) and I=I_k U I_r.}

    isSystematic, I, Tup := IsSystematic(C);
            
    return I, Tup;

end intrinsic;





/****************************************************************/
/*                                                              */
/* Function name: AllInformationSets                            */
/* Parameters: C                                                */
/* Function description: Given a systematic q-ary code C of     */
/* length n, return all the possible information sets of C as a */
/* (sorted) sequence of sequences. Each inner sequence contains */
/* a set of t coordinate positions I C {1,...,n} such that      */
/* |C|=|C_I|=q^t, where C_I={ v_I : v in C } and v_I denote     */
/* the restriction of the vector v to the coordinates in I.     */
/* The function also returns a sequence of tuples. Each         */
/* tuple < I_k, I_r > contains two sequences of integers such   */
/* that I_k is an information set for the kernel K(C)           */
/*  and I=I_k U I_r. It is not checked whether the code is      */
/* systematic or not. If the function does not succeed in       */
/* finding an information set, it returns two empty sequences.  */
/* Input parameters description:                                */
/*   - C: A q-ary code                                          */
/* Output parameters description:                               */
/*   - A sequence of all information sets for C                 */
/*   - A sequence of tuples < I_k, I_r > where I_k is an        */
/*     information set for the kernel K(C) and I=I_k U I_r is   */
/*     an information set for C.                                */
/*                                                              */
/* Signature: (<CodeFld> C) -> SeqEnum, Tup                     */
/*                                                              */
/****************************************************************/
intrinsic AllInformationSets(C::CodeFld) -> SeqEnum, Tup
{Given a systematic q-ary code C of length n, return all the possible information sets of C as a (sorted) sequence of sequences. Each inner sequence contains a set of t coordinate positions I C [1,...,n] such that |C|=|C_I|=q^t, where C_I=[ v_I : v in C ] and v_I denote the restriction of the vector v to the coordinates in I. The function also returns a sequence of tuples. Each tuple < I_k, I_r > contains two sequences of integers such that I_k is an information set for the kernel K(C) and I=I_k U I_r. It is not checked whether the code is systematic or not. If the function does not succeed in finding an information set, it returns two empty sequences.}
    
    allInfoSets := [];
    allInfoTup := [];

    if(C`IsLinear) then
        for I in AllInformationSets(C`Kernel) do
            Append(~allInfoSets, I);
            Append(~allInfoTup, <I, []>);
        end for;
        return allInfoSets, allInfoTup;
    end if;

    n := Length(C);
    t, e := ParametersInfoSet(#C, C`BaseField);
    k := Dimension(C`Kernel); 
    r := t-k;
    V := InformationSpace(C`Kernel);
    numCosetRep := #C`CosetRepresentatives;

    // #C is not a power of q 
    if e eq 1 then
        return allInfoSets, allInfoTup;
    end if;

    // #C is 1, so C is the zero code
    if t eq 0 then 
        return [], <[],[]>;
    end if;
    
    if k ne 0 then
        for IK in AllInformationSets(C`Kernel) do
            GIK := DiagonalizeKernelFromInfoSetK(C`Kernel, IK);
            newCosetRepresentatives := [v - V![v[i] : i in IK] * GIK : v in C`CosetRepresentatives];
            for IR in Subsets({1..n} diff Set(IK), r) do
                setRepresentatives := { [v[i] : i in IR] : v in newCosetRepresentatives };
                if #setRepresentatives eq numCosetRep then
                    newI := Sort(IK cat Setseq(IR));
                    if newI notin allInfoSets then
                        Append(~allInfoSets, newI);
                    end if;
                    Append(~allInfoTup, <IK, Setseq(IR)>);
                end if;
            end for; 
        end for;
        
    else
        for IR in Subsets({1..n}, t) do
            CIcodewords := { [v[i] : i in IR] : v in C`CosetRepresentatives };
            if #CIcodewords eq #C then
                newI := Sort(Setseq(IR));
                Append(~allInfoSets, newI);
                Append(~allInfoTup, <[], newI>);
            end if;
        end for;
    end if;

    return allInfoSets, allInfoTup;

end intrinsic;






/*************************************************************************************************/

intrinsic IsSystematicV6(C::CodeFld) -> BoolElt, SeqEnum, Tup
{}

    n := Length(C);
    t, e := ParametersInfoSet(#C, C`BaseField);
    k := Dimension(C`Kernel);
    r := t-k;
    G := GeneratorMatrix(C`Kernel); 

    if e ne 1 then
        return false, [], <[],[]>;
    end if;

    if t eq 0 then 
        return true, [], <[], []>;
    end if;

    for IK in Subsets({1..n}, k) do
        if(Rank(Submatrix(G, [1..k], Setseq(IK))) eq k) then
            for IR in Subsets({1..n} diff IK, r) do
                newI := IK join IR;
                if(IsInformationSetV2(C, newI)) then
                    return true, Setseq(newI), <Setseq(IK), Setseq(IR)>;
                end if;
            end for;
        end if;
    end for;

    return false, [], <[],[]>;

end intrinsic;

intrinsic IsSystematicV5(C::CodeFld) -> BoolElt, SeqEnum, Tup
{}

    n := Length(C);
    t, e := ParametersInfoSet(#C, C`BaseField);
    k := Dimension(C`Kernel);
    r := t-k;
    G := GeneratorMatrix(C`Kernel); 

    if e ne 1 then
        return false, [], <[],[]>;
    end if;

    if t eq 0 then 
        return true, [], <[], []>;
    end if;

    for IK in Subsets({1..n}, k) do
        if(Rank(Submatrix(G, [1..k], Setseq(IK))) eq k) then
            for IR in Subsets({1..n} diff IK, r) do
                newI := IK join IR;
                if(IsInformationSetV3(C, newI)) then
                    return true, Setseq(newI), <Setseq(IK), Setseq(IR)>;
                end if;
            end for;
        end if;
    end for;

    return false, [], <[],[]>;

end intrinsic;

intrinsic IsSystematicV4(C::CodeFld) -> BoolElt, SeqEnum, Tup
{}

    n := Length(C);
    t, e := ParametersInfoSet(#C, C`BaseField);
    k := Dimension(C`Kernel);
    r := t-k;
    G := GeneratorMatrix(C`Kernel); 

    if e ne 1 then
        return false, [], <[],[]>;
    end if;

    if t eq 0 then 
        return true, [], <[], []>;
    end if;

    for IK in Subsets({1..n}, k) do
        if(Rank(Submatrix(G, [1..k], Setseq(IK))) eq k) then
            for IR in Subsets({1..n} diff IK, r) do
                newI := IK join IR;
                if(IsInformationSet(C, newI)) then
                    return true, Setseq(newI), <Setseq(IK), Setseq(IR)>;
                end if;
            end for;
        end if;
    end for;

    return false, [], <[],[]>;

end intrinsic;

intrinsic IsSystematicV3(C::CodeFld) -> BoolElt, SeqEnum
{}
    
    n := Length(C);
    t, e := ParametersInfoSet(#C, C`BaseField);
      
    if e ne 1 then
        return false, [];
    end if;

    if t eq 0 then 
        return true, [];
    end if;

    for I in Subsets({1..n}, t) do
        CIcodewords := { [(k+v)[i] : i in I] : k in C`Kernel, v in C`CosetRepresentatives };
        if #CIcodewords eq #C then
            return true, Setseq(I);
        end if;
    end for;
  
    return false, [];

end intrinsic;

intrinsic IsSystematicV2(C::CodeFld) -> BoolElt, SeqEnum
{}
    
    n := Length(C);
    t, e := ParametersInfoSet(#C, C`BaseField);
      
    if e ne 1 then
        return false, [];
    end if;

    if t eq 0 then 
        return true, [];
    end if;

    for I in Subsets({1..n}, t) do
        CIcodewords := { [c[i] : i in I] : c in Set(C) }; 
        if #CIcodewords eq #C then
            return true, Setseq(I);
        end if;
    end for;
  
    return false, [];

end intrinsic;

intrinsic IsInformationSetV3(C::CodeFld, I::SeqEnum) -> BoolElt
{}

    n := Length(C);
    require (I subset [1..n]):"Argument 2 is not a subset of coordinates";
    I := Set(I);
    t, e := ParametersInfoSet(#C, C`BaseField);

    if e ne 1 then
        return false;
    elif t ne #I then
        return false;
    else
        CIcodewords := { [c[i] : i in I] : c in Set(C) }; 

        if #CIcodewords eq #C then
            return true;
        else
            return false;
        end if; 
    end if;

end intrinsic;

intrinsic IsInformationSetV2(C::CodeFld, I::SeqEnum) -> BoolElt
{}

    n := Length(C);
    require (I subset [1..n]):"Argument 2 is not a subset of coordinates";
    I := Set(I);
    t, e := ParametersInfoSet(#C, C`BaseField);

    if e ne 1 then
        return false;
    elif t ne #I then
        return false;
    else
        CIcodewords := { [(k+v)[i] : i in I] : k in C`Kernel, 
                                               v in C`CosetRepresentatives };
        
        if #CIcodewords eq #C then
            return true;
        else
            return false;
        end if; 
    end if;

end intrinsic;

intrinsic IsInformationSetV2(C::CodeFld, I::Tup) -> BoolElt
{}

    n := Length(C);
    t, _ := ParametersInfoSet(#C, C`BaseField);
    k := Dimension(C`Kernel);  
    r := t-k;
    IK := I[1];
    IR := I[2];
    
    require (#IK eq k):"Argument 2 must contain ",k," coordinates";
    require (#IR eq r):"Argument 3 must contain ",r," coordinates";
    require (IK subset [1..n]):"Argument 2 is not a subset of coordinates";
    require (IR subset [1..n]):"Argument 3 is not a subset of coordinates";
    
    if(IsInformationSet(C`Kernel, IK)) then
        if(IsInformationSet(C, IK cat IR)) then
            return true;
        end if;
    end if;
    return false;
    
end intrinsic;

intrinsic AllInformationSetsV4(C::CodeFld) -> SeqEnum
{}
    
    n := Length(C);
    t, e := ParametersInfoSet(#C, C`BaseField);
    allInfoSets := [];

    require e eq 1:"Argument 1 is not a systematic code";

    if t eq 0 then 
        return [];
    end if;

    for I in Subsets({1..n},t) do
        CIcodewords := { [c[i] : i in I] : c in Set(C) }; 
        if #CIcodewords eq #C then
            Append(~allInfoSets, Sort(Setseq(I)));
        end if;
    end for;

    require allInfoSets ne []:"Argument 1 is not a systematic code";

    return allInfoSets;

end intrinsic;

intrinsic AllInformationSetsV3(C::CodeFld) -> SeqEnum
{}
    
    n := Length(C);
    t, e := ParametersInfoSet(#C, C`BaseField);
    allInfoSets := [];

    require e eq 1:"Argument 1 is not a systematic code";

    if t eq 0 then 
        return [];
    end if;

    for I in Subsets({1..n},t) do
        CIcodewords := { [(k+v)[i] : i in I] : k in C`Kernel, 
                                               v in C`CosetRepresentatives };
        if #CIcodewords eq #C then
            Append(~allInfoSets, Sort(Setseq(I)));
        end if;
    end for;

    require allInfoSets ne []:"Argument 1 is not a systematic code";

    return allInfoSets;

end intrinsic;

intrinsic AllInformationSetsV2(C::CodeFld) -> SeqEnum
{}
    
    n := Length(C);
    t, e := ParametersInfoSet(#C, C`BaseField);
    allInfoSets := [];

    require e eq 1:"Argument 1 is not a systematic code";

    if t eq 0 then 
        return [];
    end if;

    for I in Subsets({1..n},t) do
        newI := Setseq(I);
        if(IsInformationSet(C,newI)) then
            Append(~allInfoSets, newI);
        end if;
    end for;

    require allInfoSets ne []:"Argument 1 is not a systematic code";
    
    return allInfoSets;

end intrinsic;

