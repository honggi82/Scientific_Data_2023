% Neuromag layout of whole channels
% This function returns the positions of 306 MEG channels
% Copyright by Hong Gi Yeom
% You can use this code for academic purposes

function [position]=neuromag_layout_all()

p=[
1 -73.416206 33.416687 5.000000 4.800000 0113
2 -73.416206 38.416687 5.000000 4.800000 0112
3 -67.416206 35.916687 5.000000 4.800000 0111
4 -59.602242 38.489067 5.000000 4.800000 0122
5 -59.602242 43.489067 5.000000 4.800000 0123
6 -53.602242 40.989067 5.000000 4.800000 0121
7 -68.018288 18.676970 5.000000 4.800000 0132
8 -68.018288 23.676970 5.000000 4.800000 0133
9 -62.018288 21.176970 5.000000 4.800000 0131
10 -80.582848 8.095787 5.000000 4.800000 0143
11 -80.582848 13.095787 5.000000 4.800000 0142
12 -74.582848 10.595787 5.000000 4.800000 0141
13 -56.595154 17.019251 5.000000 4.800000 0213
14 -56.595154 22.019251 5.000000 4.800000 0212
15 -50.595154 19.519251 5.000000 4.800000 0211
16 -44.599728 17.543873 5.000000 4.800000 0222
17 -44.599728 22.543873 5.000000 4.800000 0223
18 -38.599728 20.043873 5.000000 4.800000 0221
19 -47.416420 -0.216784 5.000000 4.800000 0232
20 -47.416420 4.783216 5.000000 4.800000 0233
21 -41.416420 2.283216 5.000000 4.800000 0231
22 -59.280643 -2.761772 5.000000 4.800000 0243
23 -59.280643 2.238228 5.000000 4.800000 0242
24 -53.280643 -0.261772 5.000000 4.800000 0241
25 -39.790501 47.430138 5.000000 4.800000 0313
26 -39.790501 52.430138 5.000000 4.800000 0312
27 -33.790501 49.930138 5.000000 4.800000 0311
28 -38.014336 32.768585 5.000000 4.800000 0322
29 -38.014336 37.768585 5.000000 4.800000 0323
30 -32.014336 35.268585 5.000000 4.800000 0321
31 -27.679966 28.868065 5.000000 4.800000 0333
32 -27.679966 33.868065 5.000000 4.800000 0332
33 -21.679966 31.368065 5.000000 4.800000 0331
34 -49.684467 34.878434 5.000000 4.800000 0343
35 -49.684467 40.078434 5.000000 4.800000 0342
%35 -49.684467 39.078434 5.000000 4.800000 0342
36 -43.684467 36.578434 5.000000 4.800000 0341
37 -32.997990 15.607347 5.000000 4.800000 0413
38 -32.997990 20.607347 5.000000 4.800000 0412
39 -26.997990 18.107347 5.000000 4.800000 0411
40 -21.084751 13.953575 5.000000 4.800000 0422
41 -21.084751 18.953575 5.000000 4.800000 0423
42 -15.084751 16.453575 5.000000 4.800000 0421
43 -21.930935 -0.085500 5.000000 4.800000 0432
44 -21.930935 4.914500 5.000000 4.800000 0433
45 -15.930935 2.414500 5.000000 4.800000 0431
46 -34.824663 0.362587 5.000000 4.800000 0443
47 -34.824663 5.362587 5.000000 4.800000 0442
48 -28.824663 2.862587 5.000000 4.800000 0441
49 -27.861498 55.439636 5.000000 4.800000 0513
50 -27.861498 60.439636 5.000000 4.800000 0512
51 -21.861498 57.939636 5.000000 4.800000 0511
52 -15.506709 59.619865 5.000000 4.800000 0523
53 -15.506709 64.619865 5.000000 4.800000 0522
54 -9.506709 62.119865 5.000000 4.800000 0521
55 -14.616095 49.308380 5.000000 4.800000 0532
56 -14.616095 54.308380 5.000000 4.800000 0533
57 -8.616095 51.808380 5.000000 4.800000 0531
58 -27.240477 43.863430 5.000000 4.800000 0542
59 -27.240477 48.863430 5.000000 4.800000 0543
60 -21.240477 46.363430 5.000000 4.800000 0541
61 -14.782405 38.147827 5.000000 4.800000 0613
62 -14.782405 43.147827 5.000000 4.800000 0612
63 -8.782405 40.647827 5.000000 4.800000 0611
64 -2.967276 27.260933 5.000000 4.800000 0622
65 -2.967276 32.260933 5.000000 4.800000 0623
66 3.032724 29.760933 5.000000 4.800000 0621
67 -9.094766 14.700909 5.000000 4.800000 0633
68 -9.094766 19.700909 5.000000 4.800000 0632
69 -3.094766 17.200909 5.000000 4.800000 0631
70 -15.199021 26.631405 5.000000 4.800000 0642
71 -15.199021 31.631405 5.000000 4.800000 0643
72 -9.199021 29.131405 5.000000 4.800000 0641
73 -9.246834 1.693846 5.000000 4.800000 0713
74 -9.246834 6.693846 5.000000 4.800000 0712
75 -3.246834 4.193846 5.000000 4.800000 0711
76 3.314525 1.573887 5.000000 4.800000 0723
77 3.314525 6.573887 5.000000 4.800000 0722
78 9.314525 4.873887 5.000000 4.800000 0721
79 3.387173 -10.588106 5.000000 4.800000 0733
80 3.387173 -5.588106 5.000000 4.800000 0732
81 9.387173 -8.088106 5.000000 4.800000 0731
82 -9.422897 -10.519942 5.000000 4.800000 0743
83 -9.422897 -5.519942 5.000000 4.800000 0742
84 -3.422897 -8.019942 5.000000 4.800000 0741
85 -2.962408 61.007698 5.000000 4.800000 0813
86 -2.962408 66.007698 5.000000 4.800000 0812
%86 -2.962408 65.007698 5.000000 4.800000 0812
87 3.037592 63.507698 5.000000 4.800000 0811
88 -2.965545 50.641838 5.000000 4.800000 0822
89 -2.965545 55.641838 5.000000 4.800000 0823
90 3.034455 53.141838 5.000000 4.800000 0821
91 9.504830 59.655254 5.000000 4.800000 0913
92 9.504830 64.655254 5.000000 4.800000 0912
93 15.504830 62.155254 5.000000 4.800000 0911
94 21.967310 55.408710 5.000000 4.800000 0923
95 21.967310 60.408710 5.000000 4.800000 0922
96 27.967310 57.908710 5.000000 4.800000 0921
97 21.254196 43.889683 5.000000 4.800000 0932
98 21.254196 48.889683 5.000000 4.800000 0933
99 27.254196 46.389683 5.000000 4.800000 0931
100 8.661931 49.358044 5.000000 4.800000 0942
101 8.661931 54.358044 5.000000 4.800000 0943
102 14.661931 51.858044 5.000000 4.800000 0941
103 -2.967087 39.669956 5.000000 4.800000 1013
104 -2.967087 44.669956 5.000000 4.800000 1012
105 3.032913 42.169956 5.000000 4.800000 1011
106 8.751018 38.154079 5.000000 4.800000 1023
107 8.751018 43.154079 5.000000 4.800000 1022
108 14.751018 40.654079 5.000000 4.800000 1021
109 9.123913 26.648697 5.000000 4.800000 1032
110 9.123913 31.648697 5.000000 4.800000 1033
111 15.123913 29.148697 5.000000 4.800000 1031
112 3.200539 14.795620 5.000000 4.800000 1043
113 3.200539 19.795620 5.000000 4.800000 1042
114 9.200539 17.295620 5.000000 4.800000 1041
115 15.014965 13.912239 5.000000 4.800000 1112
116 15.014965 18.912239 5.000000 4.800000 1113
117 21.014965 16.412239 5.000000 4.800000 1111
118 26.958527 15.562130 5.000000 4.800000 1123
119 26.958527 20.562130 5.000000 4.800000 1122
120 32.958527 18.062130 5.000000 4.800000 1121
121 28.757563 0.227141 5.000000 4.800000 1133
122 28.757563 5.227141 5.000000 4.800000 1132
123 34.757563 2.727141 5.000000 4.800000 1131
124 15.882982 0.037700 5.000000 4.800000 1142
125 15.882982 5.037700 5.000000 4.800000 1143
126 21.882982 2.537700 5.000000 4.800000 1141
127 33.958897 47.388790 5.000000 4.800000 1213
128 33.958897 52.388790 5.000000 4.800000 1212
129 39.958897 49.888790 5.000000 4.800000 1211
130 43.923473 33.914738 5.000000 4.800000 1223
131 43.923473 38.914738 5.000000 4.800000 1222
132 49.923473 36.414738 5.000000 4.800000 1221
133 32.014336 32.768585 5.000000 4.800000 1232
134 32.014336 37.768585 5.000000 4.800000 1233
135 38.014336 35.268585 5.000000 4.800000 1231
136 21.600079 28.898149 5.000000 4.800000 1243
137 21.600079 33.898149 5.000000 4.800000 1242
138 27.600079 31.398149 5.000000 4.800000 1241
139 38.599728 17.543867 5.000000 4.800000 1312
140 38.599728 22.543867 5.000000 4.800000 1313
141 44.599728 20.043867 5.000000 4.800000 1311
142 50.558392 16.887651 5.000000 4.800000 1323
143 50.558392 21.887651 5.000000 4.800000 1322
144 56.558392 19.387651 5.000000 4.800000 1321
145 53.420483 -2.919475 5.000000 4.800000 1333
146 53.420483 2.080525 5.000000 4.800000 1332
147 59.420483 -0.419475 5.000000 4.800000 1331
148 41.371586 -0.216817 5.000000 4.800000 1342
149 41.371586 4.783183 5.000000 4.800000 1343
150 47.371586 2.283183 5.000000 4.800000 1341
151 53.704369 38.563030 5.000000 4.800000 1412
152 53.704369 43.563030 5.000000 4.800000 1413
153 59.704369 41.063030 5.000000 4.800000 1411
154 67.119286 33.843739 5.000000 4.800000 1423
155 67.119286 38.843739 5.000000 4.800000 1422
156 73.119286 36.343739 5.000000 4.800000 1421
157 74.438919 8.335863 5.000000 4.800000 1433
158 74.438919 13.335863 5.000000 4.800000 1432
159 80.438919 10.835863 5.000000 4.800000 1431
160 61.883209 18.562304 5.000000 4.800000 1442
161 61.883209 23.562304 5.000000 4.800000 1443
162 67.883209 21.062304 5.000000 4.800000 1441
163 -71.298943 -4.707253 5.000000 4.800000 1512
164 -71.298943 0.292747 5.000000 4.800000 1513
165 -65.298943 -2.207253 5.000000 4.800000 1511
166 -67.281609 -25.407852 5.000000 4.800000 1522
167 -67.281609 -20.407852 5.000000 4.800000 1523
168 -61.281609 -22.907852 5.000000 4.800000 1521
169 -71.702820 -40.152336 5.000000 4.800000 1533
170 -71.702820 -35.152336 5.000000 4.800000 1532
171 -65.702820 -37.652336 5.000000 4.800000 1531
172 -79.907913 -17.418098 5.000000 4.800000 1543
173 -79.907913 -12.418098 5.000000 4.800000 1542
174 -73.907913 -14.918098 5.000000 4.800000 1541
175 -56.916454 -20.312164 5.000000 4.800000 1613
176 -56.916454 -15.312164 5.000000 4.800000 1612
177 -50.916454 -17.812164 5.000000 4.800000 1611
178 -45.631779 -16.320436 5.000000 4.800000 1622
179 -45.631779 -11.320436 5.000000 4.800000 1623
180 -39.631779 -13.820436 5.000000 4.800000 1621
181 -37.896103 -30.578358 5.000000 4.800000 1632
182 -37.896103 -25.578358 5.000000 4.800000 1633
183 -31.896103 -28.078358 5.000000 4.800000 1631
184 -48.859089 -36.176094 5.000000 4.800000 1643
185 -48.859089 -31.176094 5.000000 4.800000 1642
186 -42.859089 -33.676094 5.000000 4.800000 1641
187 -56.796040 -59.082275 5.000000 4.800000 1713
188 -56.796040 -53.882275 5.000000 4.800000 1712
%188 -56.796040 -54.882275 5.000000 4.800000 1712
189 -50.796040 -56.582275 5.000000 4.800000 1711
190 -57.188797 -44.857373 5.000000 4.800000 1722
191 -57.188797 -39.057373 5.000000 4.800000 1723
192 -51.188797 -41.557373 5.000000 4.800000 1721
193 -41.902962 -58.279526 5.000000 4.800000 1732
194 -41.902962 -53.279526 5.000000 4.800000 1733
195 -35.902962 -55.779526 5.000000 4.800000 1731
196 -37.408134 -72.449036 5.000000 4.800000 1743
197 -37.408134 -67.449036 5.000000 4.800000 1742
198 -31.408134 -69.949036 5.000000 4.800000 1741
199 -33.801163 -13.768716 5.000000 4.800000 1813
200 -33.801163 -8.768716 5.000000 4.800000 1812
201 -27.801163 -11.268716 5.000000 4.800000 1811
202 -21.685101 -12.619589 5.000000 4.800000 1822
203 -21.685101 -7.619589 5.000000 4.800000 1823
204 -15.685101 -10.119589 5.000000 4.800000 1821
205 -9.600111 -22.190945 5.000000 4.800000 1832
206 -9.600111 -17.190945 5.000000 4.800000 1833
207 -3.600111 -19.690945 5.000000 4.800000 1831
208 -24.483526 -26.850609 5.000000 4.800000 1843
209 -24.483526 -21.850609 5.000000 4.800000 1842
210 -18.483526 -24.350609 5.000000 4.800000 1841
211 -25.866816 -40.850040 5.000000 4.800000 1912
212 -25.866816 -35.850040 5.000000 4.800000 1913
213 -19.866816 -38.350040 5.000000 4.800000 1911
214 -20.513481 -56.355225 5.000000 4.800000 1923
215 -20.513481 -51.355225 5.000000 4.800000 1922
216 -14.513481 -53.855225 5.000000 4.800000 1921
217 -23.428471 -67.375893 5.000000 4.800000 1932
218 -23.428471 -62.375893 5.000000 4.800000 1933
219 -17.428471 -64.875893 5.000000 4.800000 1931
220 -36.237587 -48.444530 5.000000 4.800000 1943
221 -36.237587 -43.444530 5.000000 4.800000 1942
222 -30.237587 -45.944530 5.000000 4.800000 1941
223 -10.441930 -34.308243 5.000000 4.800000 2013
224 -10.441930 -29.308243 5.000000 4.800000 2012
225 -4.441930 -31.808243 5.000000 4.800000 2011
226 4.357624 -34.289736 5.000000 4.800000 2023
227 4.357624 -29.289736 5.000000 4.800000 2022
228 10.357624 -31.789736 5.000000 4.800000 2021
229 4.645295 -46.290749 5.000000 4.800000 2032
230 4.645295 -41.290749 5.000000 4.800000 2033
231 10.645295 -43.790749 5.000000 4.800000 2031
232 -10.645079 -46.244335 5.000000 4.800000 2042
233 -10.645079 -41.244335 5.000000 4.800000 2043
234 -4.645079 -43.744335 5.000000 4.800000 2041
235 -3.052351 -58.889515 5.000000 4.800000 2113
236 -3.052351 -53.889515 5.000000 4.800000 2112
237 2.947649 -56.389515 5.000000 4.800000 2111
238 -2.999999 -70.362061 5.000000 4.800000 2122
239 -2.999999 -65.362061 5.000000 4.800000 2123
240 3.000001 -67.862061 5.000000 4.800000 2121
241 8.918572 -79.441826 5.000000 4.800000 2133
242 8.918572 -74.441826 5.000000 4.800000 2132
243 14.918572 -76.941826 5.000000 4.800000 2131
244 -14.987089 -79.428932 5.000000 4.800000 2143
245 -14.987089 -74.428932 5.000000 4.800000 2142
246 -8.987089 -76.928932 5.000000 4.800000 2141
247 15.641460 -12.579389 5.000000 4.800000 2212
248 15.641460 -7.579389 5.000000 4.800000 2213
249 21.641460 -10.079389 5.000000 4.800000 2211
250 27.786499 -13.669980 5.000000 4.800000 2223
251 27.786499 -8.669980 5.000000 4.800000 2222
252 33.786499 -11.169980 5.000000 4.800000 2221
253 18.501518 -26.949615 5.000000 4.800000 2233
254 18.501518 -21.949615 5.000000 4.800000 2232
255 24.501518 -24.449615 5.000000 4.800000 2231
256 3.641699 -22.206125 5.000000 4.800000 2242
257 3.641699 -17.206125 5.000000 4.800000 2243
258 9.641699 -19.706125 5.000000 4.800000 2241
259 19.852789 -40.871220 5.000000 4.800000 2312
260 19.852789 -35.871220 5.000000 4.800000 2313
261 25.852789 -38.371220 5.000000 4.800000 2311
262 30.078903 -48.474960 5.000000 4.800000 2323
263 30.078903 -43.474960 5.000000 4.800000 2322
264 36.078903 -45.974960 5.000000 4.800000 2321
%264 35.078903 -45.974960 5.000000 4.800000 2321
265 17.363274 -67.365387 5.000000 4.800000 2332
266 17.363274 -62.365387 5.000000 4.800000 2333
267 23.363274 -64.865387 5.000000 4.800000 2331
268 14.329920 -56.380260 5.000000 4.800000 2343
269 14.329920 -51.380260 5.000000 4.800000 2342
270 20.329920 -53.880260 5.000000 4.800000 2341
271 39.644810 -16.175139 5.000000 4.800000 2412
272 39.644810 -11.175139 5.000000 4.800000 2413
273 45.644810 -13.675139 5.000000 4.800000 2411
274 50.812263 -20.401899 5.000000 4.800000 2423
275 50.812263 -15.401899 5.000000 4.800000 2422
276 56.812263 -17.901899 5.000000 4.800000 2421
277 42.694180 -36.278580 5.000000 4.800000 2433
278 42.694180 -31.278580 5.000000 4.800000 2432
279 48.694180 -33.778580 5.000000 4.800000 2431
280 31.896111 -30.578348 5.000000 4.800000 2442
281 31.896111 -25.578348 5.000000 4.800000 2443
282 37.896111 -28.078348 5.000000 4.800000 2441
283 35.812634 -58.300888 5.000000 4.800000 2512
284 35.812634 -53.300888 5.000000 4.800000 2513
285 41.812634 -55.800888 5.000000 4.800000 2511
286 51.171906 -43.981274 5.000000 4.800000 2522
287 51.171906 -38.981274 5.000000 4.800000 2523
288 57.171906 -41.481274 5.000000 4.800000 2521
289 50.704624 -59.132656 5.000000 4.800000 2533
290 50.704624 -54.132656 5.000000 4.800000 2532
291 56.704624 -56.632656 5.000000 4.800000 2531
292 31.320171 -72.484848 5.000000 4.800000 2543
293 31.320171 -67.484848 5.000000 4.800000 2542
294 37.320171 -69.984848 5.000000 4.800000 2541
295 65.137360 -4.702045 5.000000 4.800000 2612
296 65.137360 0.297955 5.000000 4.800000 2613
297 71.137360 -2.202045 5.000000 4.800000 2611
298 73.822243 -17.329140 5.000000 4.800000 2623
299 73.822243 -12.329140 5.000000 4.800000 2622
300 79.822243 -14.829140 5.000000 4.800000 2621
301 65.490112 -40.332645 5.000000 4.800000 2633
302 65.490112 -35.332645 5.000000 4.800000 2632
303 71.490112 -37.832645 5.000000 4.800000 2631
304 61.220192 -25.385981 5.000000 4.800000 2642
305 61.220192 -20.385981 5.000000 4.800000 2643
306 67.220192 -22.885981 5.000000 4.800000 2641
];

position=p(:,2:5);
for i=3:3:306
    position(i,1)=position(i,1)-1;
end
position(:,1)=(position(:,1)-min(position(:,1))+2)/(4.3*std(position(:,1)));
position(:,2)=(position(:,2)-min(position(:,2))+6)/(4.3*std(position(:,2)));
position(:,3)=0.028;
position(:,4)=0.028;