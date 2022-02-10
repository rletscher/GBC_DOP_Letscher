function [f,fx,fxx] = neglogpost_DOP(x,parm)
%
% [f,fx,fxx] = neglogpost(x,parm);
% 
% output:
% f:   -log(prob(x|data))
% fx:  1st derivative of f w.r.t. x
% fxx: 2nd derivative of f w.r.t. x 
%
% input:
% x:   model parameters
% parm: non-optimizable parameters and data that are needed to set up the model
tic
M3d  = parm.M3d;
iwet = find(M3d(:));
nwet = length(iwet);


po4obs = parm.po4obs;
DOPobs = parm.DOPobs;

% prescribe some constants.
dpa = 365;      % days per year;
spd = 24*60^2;  % seconds per day;
spa = dpa*spd;  % seconds per year;


ipo4 = find(~isnan(po4obs(iwet))); % index for valid PO4 measurements
idop = find(~isnan(DOPobs(iwet))); % index for valid DOP measurements


% extract valid data;
PO4   = po4obs(iwet(ipo4));
DOP_o = DOPobs(iwet(idop));
parm.DOPstd = std(DOP_o);


% get P model parameters;
xp = x(1:30);
ip = [1,2,3,4,5,6,7,8,9,10,...
     11,12,13,14,15,16,17,18,19,20,...
     21,22,23,24,25,26,27,28,29,30];
% solve the steady-state equilibrium P-cycle model
[P,Px,Pxx,parm] = eqPcycle_DOP(parm,ip,xp); 
parm.P = P;     % vector for P field 
parm.Px = Px;   % first derivatives;
parm.Pxx = Pxx; % second derivatives;

% reshape vector state variables into 3d fields
DIP = M3d+nan;  DIP(iwet) = P(1:nwet);
POP = M3d+nan;  POP(iwet) = P(1+nwet:2*nwet);
DOP = M3d+nan;  DOP(iwet) = P(1+2*nwet:end);
parm.DIP = DIP(iwet); parm.POP = POP(iwet); parm.DOP = DOP(iwet);


% save P fields;
fname =sprintf('output'); % output.mat holds eq. 3D DIP, DOP, POP fields
save(fname,'DIP','DOP','POP')
% covariance matrix is diagonal with variances given by
% the inverse of the grid box volume times the variance of the obs
dVt = parm.dVt;
Wp = d0(dVt(iwet(ipo4))./(parm.DIPstd.^2*sum(dVt(iwet))));
Wo = d0(dVt(iwet(idop))./(parm.DOPstd.^2*sum(dVt(iwet))));

% 
px = zeros(nwet,30);
px(:,1:30) = Px(1:nwet,:);
pos = parm.pos;
pxx = zeros(nwet,pos(end,end));


    posP = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30;...
            0 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59;...
            0 0 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81 82 83 84 85 86 87;...
            0 0 0 88 89 90 91 92 93 94 95 96 97 98 99 100 101 102 103 104 105 106 107 108 109 110 111 112 113 114;...
            0 0 0 0 115 116 117 118 119 120 121 122 123 124 125 126 127 128 129 130 131 132 133 134 135 136 137 138 139 140;...
            0 0 0 0 0 141 142 143 144 145 146 147 148 149 150 151 152 153 154 155 156 157 158 159 160 161 162 163 164 165;...
            0 0 0 0 0 0 166 167 168 169 170 171 172 173 174 175 176 177 178 179 180 181 182 183 184 185 186 187 188 189;...
            0 0 0 0 0 0 0 190 191 192 193 194 195 196 197 198 199 200 201 202 203 204 205 206 207 208 209 210 211 212;...
            0 0 0 0 0 0 0 0 213 214 215 216 217 218 219 220 221 222 223 224 225 226 227 228 229 230 231 232 233 234;...
            0 0 0 0 0 0 0 0 0 235 236 237 238 239 240 241 242 243 244 245 246 247 248 249 250 251 252 253 254 255;...
            0 0 0 0 0 0 0 0 0 0 256 257 258 259 260 261 262 263 264 265 266 267 268 269 270 271 272 273 274 275;...
            0 0 0 0 0 0 0 0 0 0 0 276 277 278 279 280 281 282 283 284 285 286 287 288 289 290 291 292 293 294;...
            0 0 0 0 0 0 0 0 0 0 0 0 295 296 297 298 299 300 301 302 303 304 305 306 307 308 309 310 311 312;...
            0 0 0 0 0 0 0 0 0 0 0 0 0 313 314 315 316 317 318 319 320 321 322 323 324 325 326 327 328 329;...
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 330 331 332 333 334 335 336 337 338 339 340 341 342 343 344 345;...
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 346 347 348 349 350 351 352 353 354 355 356 357 358 359 360;...
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 361 362 363 364 365 366 367 368 369 370 371 372 373 374;...
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 375 376 377 378 379 380 381 382 383 384 385 386 387;...
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 388 389 390 391 392 393 394 395 396 397 398 399;...
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 400 401 402 403 404 405 406 407 408 409 410;...
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 411 412 413 414 415 416 417 418 419 420;...
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 421 422 423 424 425 426 427 428 429;...
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 430 431 432 433 434 435 436 437;...
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 438 439 440 441 442 443 444;...
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 445 446 447 448 449 450;...
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 451 452 453 454 455;...
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 456 457 458 459;...
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 460 461 462;...
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 463 464;...
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 465];

% pxx is an array whose columns correspond to the 2nd derivatives
% of the P-cycle model state with respect to the 30 P-cycle model
% parameters that we want to optimize
for i1 = 1:30
    for i2 = i1:30
        pxx(:,pos(i1,i2)) = Pxx(1:nwet,posP(i1,i2));
    end
end

Ox = Px(2*nwet+1:end,:); Oxx = Pxx(2*nwet+1:end,:);
P = DIP(iwet(ipo4)); px = px(ipo4,:); pxx = pxx(ipo4,:);
O = DOP(iwet(idop)); ox = Ox(idop,:); oxx = Oxx(idop,:);

ep = P-PO4;   % DIP error beteen mod. and obs.;
eo = O-DOP_o; % DOP error beteen mod. and obs.;


f = (0.5*(ep.'*Wp*ep + eo.'*Wo*eo)); % cost function;
fprintf('cost function value f = %3.3e \n',f);

% first derivative towards parameters;
if (nargout>1)
      fx = ep.'*Wp*px + eo.'*Wo*ox;
end

% second derivatives of f with respect to the parameters 
% computed using the chain rule
if (nargout>2)
    fxx = zeros(30,30);
    for i1 = 1:30
        for i2 = i1:30
            fxx(i1,i2) = px(:,i1).'*Wp*px(:,i2)+...
                ep.'*Wp*pxx(:,pos(i1,i2)) +...
                ox(:,i1).'*Wo*ox(:,i2)+...
                eo.'*Wo*oxx(:,pos(i1,i2));
            fxx(i2,i1) = fxx(i1,i2); % symmetric
        end
    end    
end
toc
