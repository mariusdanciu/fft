
use std::{ops, fmt::Display, f32::consts::PI};

#[derive(Debug, Clone)]
struct Complex {
    re: f32,
    im: f32
}

impl Display for Complex {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if self.im < 0f32 {
            write!(f, "{} - {}i", self.re, -1f32*self.im)
        } else {
            write!(f, "{} + {}i", self.re, self.im)
        }
    }
}

impl Complex {
    fn conjugate(&self) -> Complex {
        Complex {
            re: self.re,
            im: -1f32*self.im
        }
    }

    fn abs(&self) -> f32 {
        f32::sqrt(self.re*self.re + self.im*self.im)
    }

}

impl ops::Add<Complex> for Complex {
    type Output = Complex;

    fn add(self, rhs: Complex) -> Self::Output {
        Complex{
            re: self.re + rhs.re,
            im: self.im + rhs.im
        }
    }
}

impl ops::Sub<Complex> for Complex {
    type Output = Complex;

    fn sub(self, rhs: Complex) -> Self::Output {
        Complex{
            re: self.re - rhs.re,
            im: self.im - rhs.im
        }
    }
}

impl ops::Mul<Complex> for Complex {
    type Output = Complex;

    fn mul(self, rhs: Complex) -> Self::Output {
        Complex{
            re: self.re * rhs.re - self.im*rhs.im,
            im: self.re * rhs.im + self.im*rhs.re
        }
    }
}

impl ops::Mul<f32> for Complex {
    type Output = Complex;

    fn mul(self, rhs: f32) -> Self::Output {
        Complex{
            re: self.re * rhs,
            im: self.im * rhs
        }
    }
}

fn complex(re: f32, im:f32) -> Complex {
    Complex { re, im }
}

fn fft(x: &Vec<Complex>) -> Vec<Complex> {
    let n = x.len();

    if n == 1 {
        return vec![x[0].clone()]
    }

    let mut even: Vec<Complex> = Vec::with_capacity(n/2);
    for i in 0 ..n/2 {
        even.push(x[2*i].clone());
    }

    let even_fft = fft(&even);

    let mut odd: Vec<Complex> = Vec::with_capacity(n/2);
    for i in 0 ..n/2 {
        odd.push(x[2*i + 1].clone());
    }

    let odd_fft = fft(&odd);

    let mut y = vec![complex(0f32, 0f32); n];

    for k in 0..n/2 {
        let fact = -2f32 * PI * (k as f32)/(n as f32);
        let wt = complex(f32::cos(fact), f32::sin(fact));


        y[k] = even_fft[k].clone() + wt.clone() * odd_fft[k].clone();
        y[k + n/2] = even_fft[k].clone() - wt.clone() * odd_fft[k].clone();
    }

    y
}

fn ifft(x: &Vec<Complex>) -> Vec<Complex> {

    let n = x.len();

    let conjugates: Vec<Complex> = x.iter().map(Complex::conjugate).collect();

    let fft = fft(&conjugates);

    fft.iter().map(Complex::conjugate).map(|c| c * (1f32/(n as f32))).collect()

}

fn main() {
    println!("Fast Fourier Transform");

    let sr = 64;
    let n = 64;
    let fact = 1f32 / sr as f32;

    let mut t: Vec<f32> = Vec::with_capacity(n);
    let mut v = 0f32;
    for _ in 0..n {
        t.push(v);
        v += fact;
    }

    let mut signal: Vec<f32> = Vec::with_capacity(n);

    for k in 0..n {
        signal.push(3f32*f32::sin(2f32*PI*1f32*t[k]));
    }
   

    for k in 0..n {
        signal[k] = signal[k] + f32::sin(2f32*PI*4f32*t[k]);
    }

    for k in 0..n {
        signal[k] = signal[k] + 0.5f32*f32::sin(2f32*PI*7f32*t[k]);
    }

    for c in signal.clone() {
        println!("{}", c)
    }
    let complex_signal: Vec<Complex> = signal.iter().map(|v|complex(*v, 0f32)).collect();

    
    let fft_ampl = fft(&complex_signal);

    let ampl: Vec<f32> = fft_ampl.iter().map(Complex::abs).collect();
    let mut freqs: Vec<f32> = Vec::with_capacity(n);

    let t = n / sr;

    for i in 0..n {
        freqs.push(i as f32 / t as f32);
    }

    for i in 0..n {
        println!("Freq {} Hz, Ampl {}", freqs[i], ampl[i])
    }

    let inv = ifft(&fft_ampl);

    for i in 0..n {
        println!("{},", inv[i].re)
    }
}
