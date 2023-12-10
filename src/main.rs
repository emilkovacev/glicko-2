// Inspiration for this mini-project came from https://dominion.games/
//
// The online version of the board game uses an implementation 
// of Glicko-2 to calculate its global rankings
//
// For more details, see http://www.glicko.net/glicko/glicko2.pdf

use std::f32::consts;

const TAU: f32 = 0.5;
const EPSILON: f32 = 0.000001;
const WIN: f32 = 1.0;
const DRAW: f32 = 0.5;
const LOSS: f32 = 0.0;


// Each player has 3 attributes:
//   - r: rating
//   - RD: rating deviation
//   - sigma: rating volatility
//
#[derive(Debug, Clone)]
struct Player (f32, f32, f32);

// A player who has never played defaults to certain attributes (Glickman, 2)
fn new_player() -> Player {
    Player(1500.0, 350.0, 0.06)
}

fn glicko_2(player: Player, opps: &Vec<Player>, scores: &Vec<f32>) -> Player {
    fn convert_scale(p: &Player) -> Player {
        Player((p.0 - 1500.0) / 173.7178, p.1 / 173.7178, p.2)
    }

    // g() and E() used in variance calculations
    fn g(phi: f32) -> f32 {
        1.0 / ((1.0 + (3.0 * (phi).powf(2.0)) / (consts::PI).powf(2.0)).sqrt())
    }
    fn E(mu: f32, mu_j: f32, phi_j: f32) -> f32 {
        1.0 / (1.0 + consts::E.powf(-1.0 * g(phi_j) * (mu - mu_j)))
    }

    // Used to converge on new volatility, based on "Illinois algorithm" (Glickman, 3)
    fn f(x: f32, a: f32, delta: f32, phi: f32, v: f32) -> f32 {
        let delta_sq = delta.powf(2.0);
        let phi_sq = phi.powf(2.0);
        ((consts::E.powf(x) * (delta_sq - phi_sq - v - consts::E.powf(x))) / ((2.0 * phi_sq + v + consts::E.powf(x)).powf(2.0))) - ((x - a) / (TAU.powf(2.0)))
    }

    // update from glicko-1 to glicko-2 rating scale
    let player = convert_scale(&player);
    let Player(mu, phi, sigma) = player; 
    let mut opps = opps.clone();
    for j in 0..opps.len() {
        opps[j] = convert_scale(&opps[j].clone()); // borrow checker moment
    }

    // compute estimated variance
    let mut v = 0.0;
    for j in 0..opps.len() {
        let Player(mu_j, phi_j, _sigma_j) = opps[j];
        let s_j = scores[j];
        let g_loc = g(phi_j);
        let E_loc = E(mu, mu_j, phi_j);
        v += g_loc.powf(2.0) * E_loc * (1.0 - E_loc);
    }
    v = v.powf(-1.0);

    // calculate delta
    let mut delta = 0.0;
    for j in 0..opps.len() {
        let Player(mu_j, phi_j, _sigma_j) = opps[j];
        let s_j = scores[j];
        let g_loc = g(phi_j);
        let E_loc = E(mu, mu_j, phi_j);
        delta += v * g_loc * (s_j - E_loc);
    }

    // determine new volatility
    let a = (sigma.powf(2.0)).ln();

    let mut A = a;
    let mut B: f32;
    let delta_sq = delta.powf(2.0);
    let phi_sq = phi.powf(2.0);
    if delta_sq > phi_sq + v {
        B = (delta_sq - phi_sq - v).ln();
    } else {
        let mut k = 1.0;
        let x = a - k * TAU;
        while f(x, a, delta, phi, v) < 0.0 {
            k += 1.0;
        }
        B = a - k * TAU;
    }
    let mut f_a = f(A, a, delta, phi, v);
    let mut f_b = f(B, a, delta, phi, v);

    while (B - A).abs() > EPSILON {
        let C = A + (A - B) * f_a / (f_b - f_a);
        let f_c = f(C, a, delta, phi, v);
        if f_c * f_b <= 0.0 {
            A = B;
            f_a = f_b;
        } else {
            f_a = f_a / 2.0
        }
        B = C;
        f_b = f_c;
    }
    let sigma_prime = consts::E.powf(A / 2.0);

    // update deviation
    let phi_star = (phi_sq + sigma_prime.powf(2.0)).sqrt();
    
    // update rating and RD
    let phi_prime = 1.0 / (((1.0 / phi_star.powf(2.0)) + (1.0 / v)).sqrt());
    let mu_prime = mu + phi_prime.powf(2.0) * delta / v;

    let r_prime = 173.7178 * mu_prime + 1500.0;
    let rd_prime = 173.7178 * phi_prime;

    Player(r_prime, rd_prime, sigma_prime)
}

// Example using data from sample calculation in paper (Glickman, 4)
fn main() {
    let p1 = Player(1500.0, 200.0, 0.06);
    let o1 = Player(1400.0, 30.0, 0.06);
    let o2 = Player(1550.0, 100.0, 0.06);
    let o3 = Player(1700.0, 300.0, 0.06);

    let mut opps = vec![o1, o2, o3];
    let scores = vec![WIN, LOSS, LOSS];
    let p1 = glicko_2(p1, &mut opps, &scores);
    println!("({}, {}, {})", p1.0, p1.1, p1.2);
}


