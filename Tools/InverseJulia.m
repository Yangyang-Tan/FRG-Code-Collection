(* ::Package:: *)

(* Wolfram Language Package *)

BeginPackage["InverseJulia`"]
(* Exported symbols added here with SymbolName::usage *)  
InverseJuliaForm

Begin["`Private`"] (* Begin Private Context *) 
ParenthesesQ[str_String] := 
 StringCount[str, "("] === StringCount[str, ")"]
NumericFunctionNames = 
  ToLowerCase[
   Select[Names["System`*"], 
    MemberQ[Attributes[#], NumericFunction] &]];

InverseJuliaForm[expr_, option_Symbol : Identity, 
  rule_List : {}] := 
 option[ReleaseHold[
   Rationalize@
    ToExpression[
     StringReplace[
       rule /. Rule[a_, b_] :> Rule[ToString[a], ToString[b]]]@
      FixedPoint[
       StringReplace[
        Longest[namedfunc : (LetterCharacter .. ~~ 
              WordCharacter ...) ~~ 
           Shortest["(" ~~ (arg1__ /; ParenthesesQ[arg1]) ~~ ")"]] :> 
         If[MemberQ[NumericFunctionNames, ToLowerCase[namedfunc]], 
          GeneralUtilities`ToTitleCase[namedfunc] ~~ "[" ~~ arg1 ~~ 
           "]", namedfunc ~~ "[" ~~ arg1 ~~ "]"]], 
       StringReplace[",)"->")"]@StringDelete[StringDelete[expr, {" ", "\n"}], " "]], 
     StandardForm, HoldForm]]]


InverseJuliaForm[expr_, option_Function : Identity, 
  rule_List : {}] := 
 option[ReleaseHold[
   Rationalize@
    ToExpression[
     StringReplace[
       rule /. Rule[a_, b_] :> Rule[ToString[a], ToString[b]]]@
      FixedPoint[
       StringReplace[
        Longest[namedfunc : (LetterCharacter .. ~~ 
              WordCharacter ...) ~~ 
           Shortest["(" ~~ (arg1__ /; ParenthesesQ[arg1]) ~~ ")"]] :> 
         If[MemberQ[NumericFunctionNames, ToLowerCase[namedfunc]], 
          GeneralUtilities`ToTitleCase[namedfunc] ~~ "[" ~~ arg1 ~~ 
           "]", namedfunc ~~ "[" ~~ arg1 ~~ "]"]], 
       StringDelete[StringDelete[expr, {" ", "\n"}], " "]], 
     StandardForm, HoldForm]]]



End[] (* End Private Context *)

EndPackage[]
