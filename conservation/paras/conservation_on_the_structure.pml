# To run this file from command line:
# > pymol conservation_on_the_structure.pml
# To run from pymol itself, open pymol and find conservation_on_the_structure.pml
# in the 'File'->'Run Script ...' dropdown menu.

# If you do not have pymol installed, check
# https://pymolwiki.org/index.php/Windows_Install
# or the links for other platforms at the bottom of that page.

load gnao1A.pdb, the_whole_thing
zoom complete=1
bg_color white
hide everything
select chainA, the_whole_thing and chain A and polymer
color white,  chainA
show cartoon, chainA
show spheres, chainA
set_color c0 = [1, 0.83, 0.17]
set_color c1 = [1, 0, 0]
set_color c2 = [0.73, 0, 0]
set_color c3 = [0.47, 0, 0]
set_color c4 = [0.21, 0, 0]
set_color c5 = [0.0625, 0.0625, 0.0625]
set_color c6 = [0.125, 0.125, 0.125]
set_color c7 = [0.1875, 0.1875, 0.1875]
set_color c8 = [0.25, 0.25, 0.25]
set_color c9 = [0.3125, 0.3125, 0.3125]
set_color c10 = [0.375, 0.375, 0.375]
set_color c11 = [0.4375, 0.4375, 0.4375]
set_color c12 = [0.5, 0.5, 0.5]
set_color c13 = [0.5625, 0.5625, 0.5625]
set_color c14 = [0.625, 0.625, 0.625]
set_color c15 = [0.6875, 0.6875, 0.6875]
set_color c16 = [0.75, 0.75, 0.75]
set_color c17 = [0.8125, 0.8125, 0.8125]
set_color c18 = [0.875, 0.875, 0.875]
set_color c19 = [0.9375, 0.9375, 0.9375]
set_color c20 = [1, 1, 1]
color c11, resid 160 and chain A
color c9, resid 338 and chain A
color c4, resid 247 and chain A
color c17, resid 135 and chain A
color c11, resid 65 and chain A
color c4, resid 190 and chain A
color c20, resid 112 and chain A
color c15, resid 114 and chain A
color c5, resid 347 and chain A
color c11, resid 238 and chain A
color c4, resid 46 and chain A
color c4, resid 149 and chain A
color c6, resid 49 and chain A
color c6, resid 54 and chain A
color c4, resid 176 and chain A
color c20, resid 119 and chain A
color c4, resid 271 and chain A
color c12, resid 128 and chain A
color c5, resid 325 and chain A
color c10, resid 144 and chain A
color c5, resid 317 and chain A
color c7, resid 208 and chain A
color c4, resid 273 and chain A
color c4, resid 51 and chain A
color c16, resid 66 and chain A
color c8, resid 69 and chain A
color c10, resid 150 and chain A
color c18, resid 105 and chain A
color c10, resid 308 and chain A
color c5, resid 217 and chain A
color c18, resid 280 and chain A
color c13, resid 225 and chain A
color c18, resid 142 and chain A
color c4, resid 45 and chain A
color c17, resid 194 and chain A
color c8, resid 41 and chain A
color c16, resid 164 and chain A
color c7, resid 289 and chain A
color c11, resid 257 and chain A
color c4, resid 201 and chain A
color c17, resid 123 and chain A
color c14, resid 64 and chain A
color c18, resid 106 and chain A
color c4, resid 203 and chain A
color c8, resid 226 and chain A
color c4, resid 278 and chain A
color c18, resid 121 and chain A
color c4, resid 55 and chain A
color c14, resid 110 and chain A
color c4, resid 187 and chain A
color c4, resid 326 and chain A
color c12, resid 159 and chain A
color c9, resid 303 and chain A
color c7, resid 162 and chain A
color c12, resid 72 and chain A
color c8, resid 175 and chain A
color c11, resid 301 and chain A
color c11, resid 192 and chain A
color c5, resid 140 and chain A
color c15, resid 297 and chain A
color c4, resid 267 and chain A
color c6, resid 88 and chain A
color c16, resid 154 and chain A
color c12, resid 282 and chain A
color c4, resid 331 and chain A
color c4, resid 37 and chain A
color c10, resid 333 and chain A
color c4, resid 44 and chain A
color c8, resid 199 and chain A
color c18, resid 169 and chain A
color c5, resid 284 and chain A
color c6, resid 152 and chain A
color c11, resid 233 and chain A
color c18, resid 98 and chain A
color c19, resid 59 and chain A
color c7, resid 231 and chain A
color c4, resid 56 and chain A
color c4, resid 136 and chain A
color c8, resid 61 and chain A
color c15, resid 189 and chain A
color c4, resid 264 and chain A
color c16, resid 294 and chain A
color c15, resid 310 and chain A
color c10, resid 221 and chain A
color c19, resid 99 and chain A
color c13, resid 58 and chain A
color c11, resid 178 and chain A
color c18, resid 126 and chain A
color c11, resid 103 and chain A
color c5, resid 206 and chain A
color c9, resid 223 and chain A
color c19, resid 96 and chain A
color c15, resid 101 and chain A
color c4, resid 157 and chain A
color c10, resid 287 and chain A
color c12, resid 86 and chain A
color c10, resid 210 and chain A
color c18, resid 89 and chain A
color c12, resid 323 and chain A
color c19, resid 306 and chain A
color c13, resid 70 and chain A
color c19, resid 292 and chain A
color c16, resid 262 and chain A
color c6, resid 275 and chain A
color c9, resid 321 and chain A
color c14, resid 259 and chain A
color c19, resid 167 and chain A
color c10, resid 95 and chain A
color c5, resid 197 and chain A
color c15, resid 240 and chain A
color c4, resid 336 and chain A
color c4, resid 254 and chain A
color c8, resid 73 and chain A
color c4, resid 182 and chain A
color c4, resid 184 and chain A
color c15, resid 252 and chain A
color c4, resid 269 and chain A
color c8, resid 340 and chain A
color c16, resid 299 and chain A
color c16, resid 236 and chain A
color c4, resid 131 and chain A
color c19, resid 133 and chain A
color c11, resid 85 and chain A
color c12, resid 290 and chain A
color c17, resid 314 and chain A
color c4, resid 235 and chain A
color c6, resid 147 and chain A
color c11, resid 349 and chain A
color c4, resid 260 and chain A
color c13, resid 81 and chain A
color c4, resid 212 and chain A
color c12, resid 249 and chain A
color c9, resid 138 and chain A
color c4, resid 335 and chain A
color c8, resid 214 and chain A
color c4, resid 91 and chain A
color c6, resid 312 and chain A
color c10, resid 68 and chain A
color c6, resid 173 and chain A
color c11, resid 305 and chain A
color c5, resid 244 and chain A
color c16, resid 108 and chain A
color c5, resid 228 and chain A
color c4, resid 48 and chain A
color c10, resid 171 and chain A
color c7, resid 219 and chain A
color c14, resid 276 and chain A
color c11, resid 342 and chain A
color c5, resid 250 and chain A
color c17, resid 94 and chain A
color c7, resid 328 and chain A
color c14, resid 77 and chain A
color c18, resid 125 and chain A
color c6, resid 319 and chain A
color c7, resid 84 and chain A
color c12, resid 180 and chain A
color c6, resid 344 and chain A
color c20, resid 117 and chain A
color c4, resid 205 and chain A
color c6, resid 242 and chain A
color c7, resid 177 and chain A
color c17, resid 71 and chain A
color c15, resid 158 and chain A
color c19, resid 316 and chain A
color c18, resid 120 and chain A
color c6, resid 200 and chain A
color c4, resid 185 and chain A
color c16, resid 111 and chain A
color c8, resid 255 and chain A
color c4, resid 216 and chain A
color c6, resid 279 and chain A
color c20, resid 113 and chain A
color c5, resid 300 and chain A
color c8, resid 288 and chain A
color c6, resid 42 and chain A
color c4, resid 143 and chain A
color c6, resid 62 and chain A
color c12, resid 274 and chain A
color c9, resid 141 and chain A
color c4, resid 246 and chain A
color c7, resid 168 and chain A
color c16, resid 330 and chain A
color c19, resid 97 and chain A
color c9, resid 198 and chain A
color c9, resid 74 and chain A
color c16, resid 346 and chain A
color c4, resid 87 and chain A
color c13, resid 230 and chain A
color c13, resid 295 and chain A
color c8, resid 265 and chain A
color c9, resid 272 and chain A
color c4, resid 38 and chain A
color c14, resid 163 and chain A
color c18, resid 296 and chain A
color c19, resid 302 and chain A
color c14, resid 191 and chain A
color c14, resid 90 and chain A
color c11, resid 161 and chain A
color c10, resid 193 and chain A
color c7, resid 266 and chain A
color c4, resid 204 and chain A
color c14, resid 345 and chain A
color c6, resid 148 and chain A
color c12, resid 239 and chain A
color c13, resid 124 and chain A
color c19, resid 122 and chain A
color c4, resid 76 and chain A
color c4, resid 202 and chain A
color c10, resid 79 and chain A
color c7, resid 80 and chain A
color c4, resid 339 and chain A
color c5, resid 304 and chain A
color c11, resid 137 and chain A
color c13, resid 245 and chain A
color c4, resid 232 and chain A
color c4, resid 270 and chain A
color c5, resid 151 and chain A
color c15, resid 107 and chain A
color c9, resid 334 and chain A
color c4, resid 215 and chain A
color c7, resid 256 and chain A
color c17, resid 309 and chain A
color c9, resid 153 and chain A
color c6, resid 83 and chain A
color c8, resid 227 and chain A
color c4, resid 93 and chain A
color c9, resid 283 and chain A
color c13, resid 186 and chain A
color c17, resid 281 and chain A
color c12, resid 332 and chain A
color c4, resid 52 and chain A
color c4, resid 209 and chain A
color c5, resid 327 and chain A
color c18, resid 315 and chain A
color c15, resid 75 and chain A
color c6, resid 234 and chain A
color c14, resid 118 and chain A
color c18, resid 129 and chain A
color c19, resid 166 and chain A
color c4, resid 293 and chain A
color c18, resid 261 and chain A
color c4, resid 322 and chain A
color c4, resid 53 and chain A
color c19, resid 92 and chain A
color c9, resid 291 and chain A
color c16, resid 196 and chain A
color c7, resid 263 and chain A
color c5, resid 224 and chain A
color c12, resid 139 and chain A
color c18, resid 67 and chain A
color c12, resid 337 and chain A
color c14, resid 248 and chain A
color c10, resid 104 and chain A
color c18, resid 102 and chain A
color c5, resid 222 and chain A
color c17, resid 82 and chain A
color c13, resid 145 and chain A
color c4, resid 348 and chain A
color c4, resid 237 and chain A
color c4, resid 324 and chain A
color c4, resid 170 and chain A
color c13, resid 78 and chain A
color c4, resid 132 and chain A
color c11, resid 127 and chain A
color c6, resid 253 and chain A
color c4, resid 329 and chain A
color c5, resid 156 and chain A
color c13, resid 115 and chain A
color c12, resid 318 and chain A
color c4, resid 207 and chain A
color c4, resid 251 and chain A
color c4, resid 47 and chain A
color c12, resid 181 and chain A
color c14, resid 286 and chain A
color c5, resid 183 and chain A
color c7, resid 50 and chain A
color c4, resid 229 and chain A
color c15, resid 218 and chain A
color c14, resid 307 and chain A
color c4, resid 134 and chain A
color c19, resid 109 and chain A
color c13, resid 311 and chain A
color c14, resid 258 and chain A
color c8, resid 277 and chain A
color c15, resid 313 and chain A
color c15, resid 285 and chain A
color c13, resid 100 and chain A
color c5, resid 36 and chain A
color c5, resid 220 and chain A
color c4, resid 39 and chain A
color c5, resid 213 and chain A
color c17, resid 63 and chain A
color c4, resid 179 and chain A
color c16, resid 116 and chain A
color c4, resid 155 and chain A
color c7, resid 211 and chain A
color c9, resid 320 and chain A
color c4, resid 40 and chain A
color c12, resid 188 and chain A
color c4, resid 57 and chain A
color c4, resid 146 and chain A
color c8, resid 60 and chain A
color c13, resid 241 and chain A
color c4, resid 174 and chain A
color c4, resid 35 and chain A
color c4, resid 243 and chain A
color c17, resid 298 and chain A
color c4, resid 268 and chain A
color c6, resid 43 and chain A
color c5, resid 343 and chain A
color c4, resid 341 and chain A
color c19, resid 165 and chain A
color c15, resid 172 and chain A
color c15, resid 195 and chain A
color c14, resid 130 and chain A
select heteroatoms,  hetatm and not solvent
select other_chains, not chain A
select struct_water, solvent and chain A
select metals, symbol  mg+ca+fe+zn+na+k+mn+cu+ni+cd+i
cartoon putty
show  cartoon,  other_chains
hide  spheres,   heteroatoms
show  sticks,   heteroatoms
show  spheres,  struct_water
show  spheres,  metals
color palecyan, struct_water
color lightteal, other_chains or heteroatoms
color magenta, metals
zoom  chain A
select poorly_scoring, resid 160+135+65+112+114+238+119+128+66+150+105
select poorly_scoring, poorly_scoring or resid 308+280+225+142+194+164+257+123+64+106+121
select poorly_scoring, poorly_scoring or resid 110+159+72+301+192+297+154+282+169+233+98
select poorly_scoring, poorly_scoring or resid 59+189+294+310+221+99+58+178+126+103+96
select poorly_scoring, poorly_scoring or resid 101+287+86+210+89+323+306+70+292+262+259
select poorly_scoring, poorly_scoring or resid 167+95+240+252+299+236+133+85+290+314+349
select poorly_scoring, poorly_scoring or resid 81+249+68+305+108+171+276+342+94+77+125
select poorly_scoring, poorly_scoring or resid 180+117+71+158+316+120+111+113+274+330+97
select poorly_scoring, poorly_scoring or resid 346+230+295+163+296+302+191+90+161+193+345
select poorly_scoring, poorly_scoring or resid 239+124+122+79+137+245+107+309+186+281+332
select poorly_scoring, poorly_scoring or resid 315+75+118+129+166+261+92+196+139+67+337
select poorly_scoring, poorly_scoring or resid 248+104+102+82+145+78+127+115+318+181+286
select poorly_scoring, poorly_scoring or resid 218+307+109+311+258+313+285+100+63+116+188
select poorly_scoring, poorly_scoring or resid 241+298+165+172+195+130
select poorly_scoring, poorly_scoring and chainA
deselect
