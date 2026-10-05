//This will contain the fragment shader in its entirety. should make life easier i think.
export const FRAG_SRC = `#version 300 es


precision highp float;
precision highp sampler2D;

uniform float cam_x;
uniform float cam_y;
uniform float cam_z;
uniform float M;
uniform float a;
uniform float tanHalfFov;   
uniform vec2 resolution;
uniform float DISC_IN;
const int   MAX_STEPS = 2500; 
const float H_STEP = 0.1;
const float R_ESCAPE = 71.0;
const float DISC_OUT = 20.0;
uniform sampler2D discTemp;
uniform sampler2D bbColor;
uniform float T_LUT_MIN;
uniform float T_LUT_MAX;
uniform float exposure;
uniform float T_peak;
uniform sampler2D starfield;            //I am aware this block could be cleaner but it works ok.
uniform vec2 aRange;

out vec4 fragColor;

struct state{
    vec4 x;
    vec4 p;
};

float x, y, z, r; //do not like this but i wrote all the functions like a fool.

//=============================================GENERAL HELPER FUNCTIONS====================================================

float solve_r(float x, float y, float z) {
    float a2 = a*a;
    float rho2 = x*x + y*y + z*z;
    float B = rho2 - a2;
    float S = sqrt(B*B + 4.0*a2*z*z);
    float r2 = (B >= 0.0) ? 0.5*(B + S)
                            : 2.0*a2*z*z / max(S - B, 1e-30); //looks fancy but its genuinely a quadratic
    return sqrt(max(r2, 0.0));
}

float f_(float M, float r, float a, float z) {
    float rcbd = r*r*r;    
    return (2.0*M*rcbd) / ((r*rcbd) + a*a*z*z);
}

float k_(float r, float a, float x, float y, float z, float p_t, float p_x, float p_y, float p_z) {
    float lt = 1.0;
    float lx = ((r*x) + (a*y)) / ((r*r) + (a*a));
    float ly = ((r*y) - (a*x)) / ((r*r) + (a*a));
    float lz = z / r;

    return (-lt*p_t) + (lx*p_x) + (ly*p_y) + (lz*p_z);
}

float l_x (float r, float x, float y, float a) {
    return (r*x + a*y) / (r*r + a*a);
}

float l_y (float r, float x, float y, float a) {
    return (r*y - a*x) / (r*r + a*a);
}

float l_z (float r, float z) {
    return (z / r);
}

//===================================================Redshift and Colour functions=========================================

float redshift(float r, float E_now, float Lz, float M, float a) {               
    float rt = sqrt(r);
    float rM = sqrt(M);
    float Omeg_keplr = rM / (rt*r + a*rM);
    //using modified metric components because theta is in equatorial plane.
    float g_tt = -(1.0 - 2.0*M / r);
    float g_tphi = -2.0*M*a / r;
    float g_phiphi = r*r + a*a + 2.0*M *a*a/r;
    float denom = -(g_tt + 2.0*g_tphi*Omeg_keplr + g_phiphi*Omeg_keplr*Omeg_keplr);
    float ut_em = inversesqrt(max(denom, 1.0e-6)); //guard, shouldnt be an issue but it extends to NaN issues.
    float b = Lz / E_now;
    return 1.0 / (ut_em*(1.0 - Omeg_keplr*b));
} //this redshift function is from the BL redner i had previously made, need to esnure the charts actually match but we'll feed
//here for now just to get it working.

vec3 sRGB_encode(vec3 u) {
    u = clamp(u, 0.0, 1.0);
    return mix(12.92*u, 1.055*pow(u, vec3(1.0/2.4)) - 0.055, step(0.0031308, u));
}

//===================================StarField Escape Functions======================================================

vec3 escapedir(float x, float y, float z, float p_t, float p_x, float p_y, float p_z) {
    float r = solve_r(x, y, z);
    float f = f_(M, r, a, z);
    float k = k_(r, a, x, y, z, p_t, p_x, p_y, p_z);
    vec3 l = vec3(l_x(r, x, y, a), l_y(r, x, y, a), l_z(r, z));
    
    return normalize(vec3(p_x, p_y, p_z) - f*k*l);
}

vec2 dirToUV(vec3 d) {
    float u = 0.5 + atan(d.y, d.x) / (2.0 * 3.141592654);
    float v = acos(clamp(d.z, -1.0, 1.0)) / 3.141592654;
    return vec2(u, v);
} 

// ==================================MATRIX IMPLEMENTATION===========================================================


mat4 build_metric() {
    float f  = f_(M, r, a, z);
    float lx = l_x(r, x, y, a);   
    float ly = l_y(r, x, y, a);
    float lz = l_z(r, z);

    return mat4(
        -1.0 + f,   f*lx,          f*ly,          f*lz,
         f*lx,      1.0 + f*lx*lx, f*lx*ly,       f*lx*lz,
         f*ly,      f*lx*ly,       1.0 + f*ly*ly, f*ly*lz,   // was f*lz*lz
         f*lz,      f*lx*lz,       f*ly*lz,       1.0 + f*lz*lz  // was 1.0 + lz*lz
    );
}

//===================================XDOT & PDOT FUNCTIONS============================================================

vec4 xdot(state s) {
    //etamunupnu - fklmu as reference for now
    float t = s.x.x;
    float x = s.x.y;
    float y = s.x.z;
    float z = s.x.w;

    float p_t = s.p.x;
    float p_x = s.p.y;
    float p_y = s.p.z;
    float p_z = s.p.w;

    float r = solve_r(x,y,z); //stage local, still need to call.
    float f = f_(M, r, a, z);
    float k = k_(r, a, x, y, z, p_t, p_x, p_y, p_z);

    return vec4(
    -p_t + f * k * 1.0,                        //lt is 1.0 in kerr-schild
    p_x - f*k*((r*x) + (a*y)) / ((r*r) + (a*a)),
    p_y - f*k*((r*y) - (a*x)) / ((r*r) + (a*a)),
    p_z - f*k*(z / r)
    );
}

//from here down is the r partial derivitives w.r.t spatial elements. delt = 0 (its trivial)

float rx_ (float x, float r, float y, float z, float a) {
    float holder = (2.0*r*r - (x*x + y*y + z*z) + a*a);
    return x*r / holder;
}

float ry_ (float y, float r, float x, float z, float a) {
    float holder = (2.0*r*r - (x*x + y*y + z*z) + a*a);
    return y*r / holder;
}

float rz_ (float r, float z, float a, float x, float y) {
    float holder = (2.0*r*r - (x*x + y*y + z*z) + a*a);
    return (r*r*z + a*a*z) / (r*holder);
}

//functions from here down are designed to make bug catching on pdot elements easier, truncated derivitives more or less LOL
//I for delx helpers, J for dely helpers and K for delz helpers

float x_helper (float x, float y, float z, float r, float a, float M) {
    float rx = rx_(x, r, y, z, a);
    float num = ((3.0*M*r*r*rx)*(r*r*r*r + a*a*z*z)) - ((4.0*r*r*r*rx)*(M*r*r*r));
    return 2.0 * (num) / ((r*r*r*r + a*a*z*z)*(r*r*r*r + a*a*z*z));
}

float I1_ (float x, float y, float z, float r, float a) {
    float rx = rx_(x, r, y, z, a);
    return ((((r + x*rx)*(r*r + a*a)) - ((2.0*r*rx)*(r*x + a*y))) / ((r*r + a*a)*(r*r + a*a)));
}

float I2_ (float x, float y, float z, float r, float a) {
    float rx = rx_(x, r, y, z, a);
    return ((((rx*y - a)*(r*r + a*a)) - ((r*y - a*x)*(2.0*r*rx))) / ((r*r + a*a)*(r*r + a*a)));
}

float I3_ (float x, float y, float z, float r, float a) {
    float rx = rx_(x, r, y, z, a);
    return ((-z*rx) / (r*r));
}

float y_helper (float x, float y, float z, float r, float a, float M) {
    float ry = ry_(y, r, x, z, a);
    float num = ((3.0*M*r*r*ry)*(r*r*r*r + a*a*z*z)) - ((4.0*r*r*r*ry)*(M*r*r*r));
    return 2.0 * (num) / ((r*r*r*r + a*a*z*z)*(r*r*r*r + a*a*z*z));
}

float J1_ (float x, float y, float z, float r, float a) {
    float ry = ry_(y, r, x, z, a);
    return ((((ry*x + a)*(r*r + a*a)) - ((r*x + a*y)*(2.0*r*ry))) / ((r*r + a*a)*(r*r + a*a)));
}

float J2_ (float x, float y, float z, float r, float a) {
    float ry = ry_(y, r, x, z, a);
    return ((((ry*y + r)*(r*r + a*a)) - ((2.0*r*ry)*(r*y - a*x))) / ((r*r + a*a)*(r*r + a*a)));
}

float J3_ (float x, float y, float z, float r, float a) {
    float ry = ry_(y, r, x, z, a);
    return ((-z*ry) / (r*r));
}

float z_helper (float x, float y, float z, float r, float a, float M) {
    float rz = rz_(r, z, a, x, y);
    float num = ((3.0*M*r*r*rz)*(r*r*r*r+ a*a*z*z)) - ((4.0*r*r*r*rz + 2.0*a*a*z)*(M*r*r*r));
    return 2.0 * (num) / (((r*r*r*r+ a*a*z*z)*(r*r*r*r+ a*a*z*z)));
}

float K1_ (float x, float y, float z, float r, float a) {
    float rz = rz_(r, z, a, x, y);
    return ((((rz*x)*(r*r + a*a) - (2.0*r*rz)*(r*x + a*y))) / ((r*r + a*a)*(r*r + a*a)));
}

float K2_ (float x, float y, float z, float r, float a) {
    float rz = rz_(r, z, a, x, y);
    return ((((rz*y)*(r*r + a*a) - (2.0*r*rz)*(r*y - a*x))) / ((r*r + a*a)*(r*r + a*a)));
}

float K3_ (float x, float y, float z, float r, float a) {
    float rz = rz_(r, z, a, x, y);
    return ((r - rz*z) / (r*r));
}

vec4 pdot (state s) {
//this one is a doozy, 0.5k^2del_if + fk(del_alphalalpha)p_alpha for reference
    float t = s.x.x;
    float x = s.x.y;
    float y = s.x.z;
    float z = s.x.w;

    float p_t = s.p.x;
    float p_x = s.p.y;
    float p_y = s.p.z;
    float p_z = s.p.w;
    
    float r  = solve_r(x,y,z); // same situation as in xdot.
    float f = f_(M, r, a, z);
    float k = k_(r, a, x, y, z, p_t, p_x, p_y, p_z);

    return vec4(
        0.0, //nice killing field symmetry here, the phi one is a little more complicated in kerr-schild
        (0.5*k*k*x_helper(x, y, z, r, a, M)) + (f * k * ((I1_(x, y, z, r, a)*p_x) + (I2_(x, y, z, r, a)*p_y) + (I3_(x, y, z, r, a)*p_z))),
        (0.5*k*k*y_helper(x, y, z, r, a, M)) + (f * k * ((J1_(x, y, z, r, a)*p_x) + (J2_(x, y, z, r, a)*p_y) + (J3_(x, y, z, r, a)*p_z))),
        (0.5*k*k*z_helper(x, y, z, r, a, M)) + (f * k * ((K1_(x, y, z, r, a)*p_x) + (K2_(x, y, z, r, a)*p_y) + (K3_(x, y, z, r, a)*p_z)))
    );    
}

struct dots {
    vec4 dx;
    vec4 dp;
};

dots dot_functions(state s) {
    dots d;
    d.dx = xdot(s);
    d.dp = pdot(s);
    return d;
}

state rk4step(state s, float h) {
    dots k1 = dot_functions(s);

    state s2;
    s2.x = s.x + 0.5*h*k1.dx;
    s2.p = s.p + 0.5*h*k1.dp;
    dots k2 = dot_functions(s2);

    state s3;
    s3.x = s.x + 0.5*h*k2.dx;
    s3.p = s.p + 0.5*h*k2.dp;
    dots k3 = dot_functions(s3);

    state s4;
    s4.x = s.x + h*k3.dx;
    s4.p = s.p + h*k3.dp;
    dots k4 = dot_functions(s4);

    s.x += (h/6.0)*(k1.dx + 2.0*k2.dx + 2.0*k3.dx + k4.dx);
    s.p += (h/6.0)*(k1.dp + 2.0*k2.dp + 2.0*k3.dp + k4.dp);

    return s;
}

//=========================TETRAD CONSTRUCTION & GRAM-SCHMIDT================================================================//

float A_hlp (float r, float a) {
    return r*r + a*a;
}

float p_sqr_hlp (float x, float y) {
    return x*x + y*y;
}

float Omega (float a) {
    float f = f_(M, r, a, z);
    float A = A_hlp(r, a);
    float p = p_sqr_hlp(x, y);
    
    return (f*a*A) / ((A*A) + f*a*a*p);
}

// these D and N factors were simple to factor for the first 2 vectors but it just became a headache later on
// so they wont get too much reuse

float D_ (float a) {
    float f = f_(M, r, a, z);
    float A = A_hlp(r, a);
    float p = p_sqr_hlp(x, y);

    return ((A*A) + f*a*a*p);
}

float N_ (float a) {
    float f = f_(M, r, a, z);
    float A = A_hlp(r, a);
    float p = p_sqr_hlp(x, y);

    return (((1.0 - f)*A*A) + f*a*a*p); 
}

vec4 ehat0_ () {
    float D = D_(a);
    float N = N_(a);
    float omeg = Omega(a);
    float mlt = sqrt(D/N);

    return vec4 (
        mlt,
        mlt*(-y*omeg),
        mlt*(x*omeg),
        0.0
    );
}
// 0 and 1 i did cleanly on paper, the others im gonna rely on computation a bit more for because theyre taing a couple pages up.
vec4 ehat1prime (float r, float x) {
    float f = f_(M, r, a, z);
    float A = A_hlp(r, a);
    float N = N_(a);
    float c1 = ((f*r*A*x) / N);
    float omeg = Omega(a);

    return vec4 (
        c1,
        1.0 - c1*omeg*y,
        c1*omeg*x,
        0.0
    );
}

vec4 ehat1 () {
    float f = f_(M, r, a, z);
    float A = A_hlp(r, a);
    float N = N_(a);
    float D = D_(a);
    float lx = l_x(r, x, y, a);
    vec4 ehp = ehat1prime(r, x);

    float dscrm = (1.0 + f*lx*lx) + ((f*f*r*r*A*A*x*x) / (D*N));
    return ehp / sqrt(dscrm);
}

vec4 ehat2prime () {
    float f = f_(M, r, a, z);
    float A = A_hlp(r, a);
    float N = N_(a);
    float D = D_(a);
    float lx = l_x(r, x, y, a);
    float omeg = Omega(a);
    vec4 ehp  = ehat1prime(r, x);
    mat4 g = build_metric();

    float dscrm = (1.0 + f*lx*lx) + ((f*f*r*r*A*A*x*x) / (D*N));

    vec4 g_ymu = g[2];                        
    vec4 v0 = vec4(1.0, -y*omeg, x*omeg, 0.0);  

    float s1 = dot(g_ymu, v0);        
    float s2 = s1 / (-(N / D));       
    vec4 s3 = s2 * v0;

    float s4 = dot(g_ymu, ehp);       
    float s5 = s4 / dscrm;
    vec4 s6 = s5 * ehp;

    vec4 dely = vec4(0.0, 0.0, 1.0, 0.0);
    return dely - s3 - s6;
}

vec4 ehat2 () {
    float f = f_(M, r, a, z);
    float A = A_hlp(r, a);
    float N = N_(a);
    float D = D_(a);
    float lx = l_x(r, x, y, a);
    float omeg = Omega(a);
    vec4 ehp = ehat1prime(r, x);
    mat4 g = build_metric();

    float dscrm = (1.0 + f*lx*lx) + ((f*f*r*r*A*A*x*x) / (D*N));

    vec4 g_ymu = g[2];
    vec4 v0 = vec4(1.0, -y*omeg, x*omeg, 0.0);
    vec4 dely = vec4(0.0, 0.0, 1.0, 0.0);

    float s1 = dot(g_ymu, v0);        
    float s2 = s1 / (-(N / D));
    float s4 = dot(g_ymu, ehp);       
    float s5 = s4 / dscrm;

    float dscrm2 = dot(g_ymu, dely) - s1*s2 - s4*s5;

    return ehat2prime() / sqrt(dscrm2);
}

vec4 ehat3prime() {
    float f = f_(M, r, a, z);
    float A = A_hlp(r, a);
    float N = N_(a);
    float D = D_(a);
    float lx = l_x(r, x, y, a);
    float omeg = Omega(a);
    vec4 ehp1 = ehat1prime(r, x);
    vec4 ehp2 = ehat2prime();
    vec4 dely = vec4(0.0, 0.0, 1.0, 0.0);
    mat4 g = build_metric();
    vec4 g_ymu = g[2];
    vec4 v0 = vec4(1.0, -y*omeg, x*omeg, 0.0);

    
    float dscrm = (1.0 + f*lx*lx) + ((f*f*r*r*A*A*x*x) / (D*N));
    float s1 = dot(g_ymu, v0);        
    float s2 = s1 / (-(N / D));
    float s4 = dot(g_ymu, ehp1);       
    float s5 = s4 / dscrm;  //A lot of abstraction makes this whole process butter
    float dscrm2 = dot(g_ymu, dely) - s1*s2 - s4*s5;
    
    vec4 g_zmu = g[3];
    vec4 delz = vec4(0.0, 0.0, 0.0, 1.0);

    float p1 = dot(g_zmu, v0);
    float p2 = (D / N)*p1;
    vec4 p3 = p2 * v0;

    float p4 = dot(g_zmu, ehp1);
    float p5 = p4 / dscrm;
    vec4 p6 = p5 * ehp1;

    float p7 = dot(g_zmu, ehp2);
    float p8 = p7 / dscrm2;
    vec4 p9 = p8*ehp2;

    return delz + p3 - p6 - p9;
}

vec4 ehat3 () {
    float f = f_(M, r, a, z);
    float A = A_hlp(r, a);
    float N = N_(a);
    float D = D_(a);
    float lx = l_x(r, x, y, a);
    float omeg = Omega(a);
    vec4 ehp1 = ehat1prime(r, x);
    vec4 ehp2 = ehat2prime();
    mat4 g = build_metric();

    float dscrm = (1.0 + f*lx*lx) + ((f*f*r*r*A*A*x*x) / (D*N));

    vec4 g_ymu = g[2];
    vec4 g_zmu = g[3];
    vec4 v0 = vec4(1.0, -y*omeg, x*omeg, 0.0);
    vec4 dely = vec4(0.0, 0.0, 1.0, 0.0);
    vec4 delz = vec4(0.0, 0.0, 0.0, 1.0);

    float s1 = dot(g_ymu, v0);        
    float s2 = s1 / (-(N / D));
    float s4 = dot(g_ymu, ehp1);       
    float s5 = s4 / dscrm;

    float dscrm2 = dot(g_ymu, dely) - s1*s2 - s4*s5;

    float p1 = dot(g_zmu, v0);
    float p2 = (D / N)*p1;
    float p4 = dot(g_zmu, ehp1);
    float p5 = p4 / dscrm;
    float p7 = dot(g_zmu, ehp2);
    float p8 = p7 / dscrm2;

    float dscrm3 = dot(g_zmu, delz) + (D/N)*p1*p1 - p4*p5 - p7*p8;
    return ehat3prime() / sqrt(dscrm3);
}

state init_momentum() {
    x = cam_x; y = cam_y; z = cam_z;
    r = solve_r(x, y, z);
    float f = f_(M, r, a, z);
    float lx = l_x(r, x, y, a);
    float ly = l_y(r, x, y, a);
    float lz = l_z(r, z);
    vec4 l = vec4(1.0, lx, ly, lz);

    mat4 tetrad = mat4(ehat0_(), ehat1(), ehat2(), ehat3());

    vec3 fwd = normalize(-vec3(x, y, z)); // camera positioning
    vec3 right = normalize(cross(vec3(0.0, 0.0, 1.0), fwd));
    vec3 up = cross(fwd, right);

    vec2 uv = (gl_FragCoord.xy / resolution) * 2.0 - 1.0;
    uv.x *= resolution.x / resolution.y;
    vec3 direction = normalize(vec3(uv.x*tanHalfFov, uv.y*tanHalfFov, -1.0)); //pinhole cam
    vec3 n_hat = direction.x*right + direction.y*up - direction.z*fwd;

    vec4 pmu = tetrad * vec4(1.0, n_hat);
    float s = dot(l, pmu);
    vec4 p_mu = vec4(-pmu.x, pmu.y, pmu.z, pmu.w) + f*s*l;

    state st;
    st.x = vec4(0.0, x, y, z);
    st.p = p_mu;
    return st;
} //eventually this can be moved cpu side but i dont want to do that rn

//==============================================================MAIN BRANCH=================================================

void main() {
    float r_plus = M + sqrt(M*M - a*a); // outer horizon
    vec3 color = vec3(0.0);
    bool done = false;
    state s = init_momentum();
    float r = solve_r(s.x.y, s.x.z, s.x.w);


    for (int i = 0; i < MAX_STEPS && !done; i++) {
        
        //ADDED, poor mans adaptive stepping.
        float h;
        if (r > 5.0) h = 3.0;
        else h = H_STEP;
        
        vec4  prev_x = s.x;
        vec4  prev_p = s.p;
        float z_prev = s.x.w;

        
        s = rk4step(s, h);
        x = s.x.y;  y = s.x.z;  z = s.x.w;
        float p_t = s.p.x, p_x = s.p.y, p_y = s.p.z, p_z = s.p.w;
        r = solve_r(x, y, z);


        //nan guard, realistically theyve hit the horizon if they land here.
        if (isnan(r) || isinf(r)) {
            color = vec3(0.0);
            done = true;
        }
        //horizon
        else if (r < r_plus) { 
            color = vec3(0.0);
            done = true;
        } 

        // Disc
        else if (z_prev * z < 0.0) {
            
            //MOVED, frac, xc, rho2, r_cross into this elif branch to save a little computation.
            float frac = z_prev / (z_prev - z);
            vec4  xc = mix(prev_x, s.x, frac);
            float rho2 = xc.y*xc.y + xc.z*xc.z;
            float r_cross = sqrt(max(rho2 - a*a, 0.0)); //interpolate crossing radius

            if (r_cross > DISC_IN && r_cross < DISC_OUT) {
                vec4 pc = mix(prev_p, s.p, frac);
                float E = -pc.x;
                float Lz = xc.y*pc.z - xc.z*pc.y;

                float u_r = (r_cross - DISC_IN) / (DISC_OUT - DISC_IN);
                float u_a = clamp((a - aRange.x) / (aRange.y - aRange.x), 0.0, 1.0);
                float T = texture(discTemp, vec2(u_r, u_a)).r * T_peak;

                float g = redshift(r_cross, E, Lz, M, a);
                float T_obs = g * T;
                float u_T = clamp((T_obs - T_LUT_MIN) / (T_LUT_MAX - T_LUT_MIN), 0.0, 1.0);
                vec3 hue = texture(bbColor, vec2(u_T, 0.5)).rgb;
                float brt = pow(T_obs / T_peak, 4.0);
                color = hue * brt * exposure;
                done = true;
            }
        }
        // 3) escape to background
        else if (r > R_ESCAPE) {
            //float x, float y, float z, float p_t, float p_x, float p_y, float p_z
            vec3 dir = escapedir(x, y, z, p_t, p_x, p_y, p_z);  
            color = texture(starfield, dirToUV(dir)).rgb * 250.0;
            done = true;
        }
    }
      
    //  ran out of steps without resolving 
    if (!done) {
        color = vec3(0.0, 1.0, 0.0);   // bright green = "MAX_STEPS hit", debug only
    }
    // luminance Reinhard (preserves hue) -> sRGB
    float Lum = dot(color, vec3(0.2126, 0.7152, 0.0722));
    if (Lum > 0.0) color *= (Lum / (1.0 + Lum)) / Lum;
    color = sRGB_encode(color);
    fragColor = vec4(color, 1.0);
    }

`