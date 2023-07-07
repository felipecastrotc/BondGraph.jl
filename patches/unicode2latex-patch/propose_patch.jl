# Reference
# http://milde.users.sourceforge.net/LUCR/Math/unimathsymbols.html

# The julia shorcuts are slightly different in some cases than the above
# https://docs.julialang.org/en/v1/manual/unicode-input/


patch = Dict{Char,String}(
    # https://en.wikipedia.org/wiki/Superscripts_and_Subscripts_(Unicode_block)
    # Numbers subscript
    '₀' => raw"_0",
    '₁' => raw"_1",
    '₂' => raw"_2",
    '₃' => raw"_3",
    '₄' => raw"_4",
    '₅' => raw"_5",
    '₆' => raw"_6",
    '₇' => raw"_7",
    '₈' => raw"_8",
    '₉' => raw"_9",
    '₊' => raw"_+",
    '₋' => raw"_-",
    '₌' => raw"_=",
    '₍' => raw"_(",
    '₎' => raw"_)",
    # Latin minuscule subscript
    'ₐ' => raw"_a",
    'ₑ' => raw"_e",
    'ₕ' => raw"_h",
    'ᵢ' => raw"_i",
    'ⱼ' => raw"_j",
    'ₖ' => raw"_k",
    'ₗ' => raw"_l",
    'ₘ' => raw"_m",
    'ₙ' => raw"_n",
    'ₒ' => raw"_o",
    'ₚ' => raw"_p",
    'ᵣ' => raw"_r",
    'ₛ' => raw"_s",
    'ₜ' => raw"_t",
    'ᵤ' => raw"_u",
    'ᵥ' => raw"_v",
    'ₓ' => raw"_x",
    # Greek minuscule subscript
    'ᵦ' => raw"_\beta",
    'ᵧ' => raw"_\gamma",
    'ᵨ' => raw"_\rho",
    'ᵩ' => raw"_\phi",
    'ᵪ' => raw"_\chi",
    # IPA subscript
    'ₔ' => raw"_{\mbox{\textschwa}}", # requires \usepackage{tipa}
    # Numbers superscript
    '⁰' => raw"^0",
    '¹' => raw"^1",
    '²' => raw"^2",
    '³' => raw"^3",
    '⁴' => raw"^4",
    '⁵' => raw"^5",
    '⁶' => raw"^6",
    '⁷' => raw"^7",
    '⁸' => raw"^8",
    '⁹' => raw"^9",
    '⁺' => raw"^+",
    '⁻' => raw"^-",
    '⁼' => raw"^=",
    '⁽' => raw"^(",
    '⁾' => raw"^)",
    'ꜝ' => raw"^!",
    'ꜛ' => raw"^{\uparrow}",
    'ꜜ' => raw"^{\downarrow}",
    # Latin minuscule superscript
    'ᵃ' => raw"^a",
    'ᵇ' => raw"^b",
    'ᶜ' => raw"^c",
    'ᵈ' => raw"^d",
    'ᵉ' => raw"^e",
    'ᶠ' => raw"^f",
    'ᵍ' => raw"^g",
    'ʰ' => raw"^h",
    'ⁱ' => raw"^i",
    'ʲ' => raw"^j",
    'ᵏ' => raw"^k",
    'ˡ' => raw"^l",
    'ᵐ' => raw"^m",
    'ⁿ' => raw"^n",
    'ᵒ' => raw"^o",
    'ᵖ' => raw"^p",
    'ʳ' => raw"^r",
    'ˢ' => raw"^s",
    'ᵗ' => raw"^t",
    'ᵘ' => raw"^u",
    'ᵛ' => raw"^v",
    'ʷ' => raw"^w",
    'ˣ' => raw"^x",
    'ʸ' => raw"^y",
    'ᶻ' => raw"^z",
    # Greek minuscule superscript
    'ᵝ' => raw"^\beta",
    'ᵞ' => raw"^\gamma",
    'ᵟ' => raw"^\delta",
    'ᵋ' => raw"^\epsilon",
    'ᶿ' => raw"^\theta",
    'ᶥ' => raw"^\iota",
    'ᶹ' => raw"^\upsilon",
    'ᵠ' => raw"^\phi",
    'ᵡ' => raw"^\chi",
    # Latin capital superscript
    'ᴬ' => raw"^A",
    'ᴮ' => raw"^B",
    'ꟲ' => raw"^C",
    'ᴰ' => raw"^D",
    'ᴱ' => raw"^E",
    'ꟳ' => raw"^F",
    'ᴳ' => raw"^G",
    'ᴴ' => raw"^H",
    'ᴵ' => raw"^I",
    'ᶦ' => raw"^I",
    'ᴶ' => raw"^J",
    'ᴷ' => raw"^K",
    'ᴸ' => raw"^L",
    'ᶫ' => raw"^L",
    'ᴹ' => raw"^M",
    'ᴺ' => raw"^N",
    'ᶰ' => raw"^N",
    'ᴼ' => raw"^O",
    'ᴾ' => raw"^P",
    'ꟴ' => raw"^Q",
    'ᴿ' => raw"^R",
    'ᵀ' => raw"^T",
    'ᵁ' => raw"^U",
    'ᶸ' => raw"^U",
    'ⱽ' => raw"^V",
    'ᵂ' => raw"^W",
    # Cyrillic superscript
    'ꚜ' => raw"^?",
    'ᵸ' => raw"^?",
    # IPA superscript
    'ᵅ' => raw"^?",
    'ᶞ' => raw"^{\mbox{\textipa{D}}}", # requires \usepackage{tipa}
    'ᵊ' => raw"^{\mbox{\textschwa}}", # requires \usepackage{tipa}
    'ᶪ' => raw"^?",
    'ᶴ' => raw"^{\mbox{\textipa{S}}}", # requires \usepackage{tipa}
    'ᶵ' => raw"^?",
    'ꭩ' => raw"^?",
    'ˀ' => raw"^{\mbox{\textipa{P}}}", # requires \usepackage{tipa}
    # Other superscript and subscript characters
    # Latin-1 Supplement
    'ª' => raw"^?",
    'º' => raw"^?",
    # Latin Extended-D
    'ꝰ' => raw"^?",
    'ꟸ' => raw"^?",
    'ꟹ' => raw"^?",
    # Latin Extended-E
    'ꭜ' => raw"^?",
    'ꭝ' => raw"^?",
    'ꭞ' => raw"^?",
    'ꭟ' => raw"^?",
    # Spacing Modifier Letters
    'ʱ' => raw"^?",
    'ʴ' => raw"^?",
    'ʵ' => raw"^?",
    'ʶ' => raw"^?",
    'ˁ' => raw"^?",
    'ˠ' => raw"^?",
    'ˤ' => raw"^?",
    # Phonetic Extensions
    'ᴭ' => raw"^?",
    'ᴯ' => raw"^?",
    'ᴲ' => raw"^?",
    'ᴻ' => raw"^?",
    'ᴽ' => raw"^?",
    'ᵄ' => raw"^?",
    'ᵆ' => raw"^?",
    'ᵌ' => raw"^?",
    'ᵑ' => raw"^?",
    'ᵓ' => raw"^?",
    'ᵚ' => raw"^?",
    'ᵸ' => raw"^?",
    'ᵎ' => raw"^?",
    'ᵔ' => raw"^?",
    'ᵕ' => raw"^?",
    'ᵙ' => raw"^?",
    'ᵜ' => raw"^?",
    # Phonetic Extensions Supplement
    'ᶛ' => raw"^?",
    'ᶝ' => raw"^?",
    'ᶟ' => raw"^?",
    'ᶡ' => raw"^?",
    'ᶢ' => raw"^?",
    'ᶣ' => raw"^?",
    'ᶤ' => raw"^?",
    'ᶧ' => raw"^?",
    'ᶨ' => raw"^?",
    'ᶩ' => raw"^?",
    'ᶬ' => raw"^?",
    'ᶭ' => raw"^?",
    'ᶮ' => raw"^?",
    'ᶯ' => raw"^?",
    'ᶱ' => raw"^?",
    'ᶲ' => raw"^?",
    'ᶳ' => raw"^?",
    'ᶶ' => raw"^?",
    'ᶷ' => raw"^?",
    'ᶺ' => raw"^?",
    'ᶼ' => raw"^?",
    'ᶽ' => raw"^?",
    'ᶾ' => raw"^?",
    # Cyrillic Extended-B
    'ꚝ' => raw"^?",
    # Georgian
    'ჼ' => raw"^?",
)


# Run unicodedict without the constant and then
merge!(unicodedict, patch)

str = unicodedict |> string;

io = open("unicodedict.txt", "w");
write(io, str);
close(io);

# Remember to put the raw"










ₔ = 1

'ₔ'
<<<<<<< HEAD













=======
>>>>>>> bg-doc/main
