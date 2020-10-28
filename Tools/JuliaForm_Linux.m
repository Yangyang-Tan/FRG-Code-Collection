(* ::Package:: *)

BeginPackage["JuliaForm`"]

JuliaForm

JuliaFormV

exp
Begin["`Private`"]

SetAttributes[JuliaForm, HoldFirst]

BracketQ[str_String] := StringCount[str, "["] === StringCount[str, "]"]

copyUnicode[expr_] := 
  Run["clip <", 
   Export["$Clipboard.temp", expr, "Text", 
    CharacterEncoding -> "Unicode"]]




JuliaForm[expr_, option_ : Identity] := 
 ToString[StringDelete[
   FixedPoint[
    StringReplace[
     Longest[namedfunc : ( LetterCharacter .. ~~ WordCharacter ...) ~~
         Shortest["[" ~~ (arg1__ /; BracketQ[arg1]) ~~ "]"]] :> 
      If[MemberQ[Attributes[namedfunc], NumericFunction], 
       ToLowerCase[namedfunc] ~~ "(" ~~ arg1 ~~ ")", 
       namedfunc ~~ "(" ~~ arg1 ~~ ")"]], 
    StringDelete[
     StringDelete[
      ToString[InputForm[option[expr] //. Exp[x_] :> exp[x]]], {" ", 
       "\n"}], " "]], {" ", "\n"}],CharacterEncoding -> "Unicode"]


JuliaFormV[expr_, option_Symbol : Identity, 
  vecfun_List : {}] := 
 copyUnicode[
  StringReplace[{"+"->" .+","-"->" .-", "*" -> " .*", "^" -> " .^", "/" -> " ./","Pi"->"pi"}]@
   StringDelete[
    StringReplace[Rule[ToString[#], ToString[#] <> "."] & /@ vecfun]@
     FixedPoint[
      StringReplace[
       Longest[namedfunc : ( 
            LetterCharacter .. ~~ WordCharacter ...) ~~ 
          Shortest["[" ~~ (arg1__ /; BracketQ[arg1]) ~~ "]"]] :> 
        If[MemberQ[Attributes[namedfunc], NumericFunction] || 
          namedfunc === "exp", 
         ToLowerCase[namedfunc] ~~ ".(" ~~ arg1 ~~ ")", 
         namedfunc ~~ "(" ~~ arg1 ~~ ")"]], 
      StringDelete[
       StringDelete[
        ToString[InputForm[option[expr] //. Exp[x_] :> exp[x]]], {" ",
          "\n"}], " "]], {" ", "\n"}]]


JuliaFormV[expr_, option_Function : Identity, 
  vecfun_List : {}] := 
 copyUnicode[
  StringReplace[{"+"->" .+","-"->" .-", "*" -> " .*", "^" -> " .^", "/" -> " ./","Pi"->"pi"}]@
   StringDelete[
    StringReplace[Rule[ToString[#], ToString[#] <> "."] & /@ vecfun]@
     FixedPoint[
      StringReplace[
       Longest[namedfunc : ( 
            LetterCharacter .. ~~ WordCharacter ...) ~~ 
          Shortest["[" ~~ (arg1__ /; BracketQ[arg1]) ~~ "]"]] :> 
        If[MemberQ[Attributes[namedfunc], NumericFunction] || 
          namedfunc === "exp", 
         ToLowerCase[namedfunc] ~~ ".(" ~~ arg1 ~~ ")", 
         namedfunc ~~ "(" ~~ arg1 ~~ ")"]], 
      StringDelete[
       StringDelete[
        ToString[InputForm[option[expr] //. Exp[x_] :> exp[x]]], {" ",
          "\n"}], " "]], {" ", "\n"}]]



End[]

EndPackage[]
