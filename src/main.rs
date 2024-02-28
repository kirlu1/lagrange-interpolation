use modinverse::modinverse;

fn main() {
    let points = vec![(0,1),(1,1),(4,1),(2,3),(3,3)];
    let lagrange_result = lagrange_interpolation(points);
    lagrange_result.debug();
}

type IntType = i32;

const CHARACTERISTIC : IntType = 5;
const POL_SIZE : usize = 256; 

type Pol = [IntType;POL_SIZE];

#[derive(Clone, Copy)]
struct Polynomial(Pol);


impl Polynomial {
    fn new() -> Self {
        Polynomial([0; POL_SIZE])
    }

    fn scalar_mul(mut self, factor : IntType) -> Self {
        for coefficient in self.iter_mut() {
            *coefficient*=factor;
        }
        self
    }

    fn debug(&self) {
        let mut output = String::new();
        let it = self.iter().enumerate();

        for (idx,c) in it.rev() {
            let addend = match (idx,c) {
                (_,0) => continue,
                (0,c) => format!("{}+",c),
                (1,1) => format!("x+"),
                (1,c) => format!("{}x+",c),
                (p,1) => format!("x^{}+",p),
                (p,c) => format!("{}x^{}+",c,p),
            };
            output.push_str(&addend)
        }
        if let None = output.pop() {
           println!("0");
        } else {
            println!("{}",output);
        }
    }
}

fn lagrange_interpolation(points : Vec<(IntType,IntType)>) -> Polynomial {
    let mut result = Polynomial::new();

    let (domain, image) : (Vec<_>, Vec<_>) = points.into_iter().unzip();

    for (i, y_i) in image.into_iter().enumerate() {
        let basis_pol = basis_polynomial(i, &domain).scalar_mul(y_i);
        //println!("{:?}", basis_pol.0);
        result=result+basis_pol;
    }

    result
}

fn basis_polynomial(i : usize, domain : &Vec<IntType>) -> Polynomial {
    let mut polynomial = Polynomial::new();
    polynomial[0]=1;

    let mut degree_one = Polynomial::new();
    degree_one[1]=1;

    domain.iter().enumerate()
    // skip index i 
    .filter(|(j, _)| i!=*j)
    .fold(polynomial, |acc, (_,val)| {

        let divisor = (domain[i]-val).rem_euclid(CHARACTERISTIC);
        let mut dividend = degree_one.clone();
        dividend[0]+=(-val).rem_euclid(CHARACTERISTIC);

        let Some(modinv) = modinverse(divisor, CHARACTERISTIC) else {
            println!("Ooops: {}",divisor);
            std::process::exit(0);
        };

        let factor= dividend.scalar_mul(modinv);

        acc*factor
    })
}

impl std::ops::Deref for Polynomial {
    type Target = Pol;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl std::ops::DerefMut for Polynomial {

    fn deref_mut(&mut self) -> &mut Pol {
        &mut self.0
    }
}

impl std::ops::Add for Polynomial {
    type Output = Self;

    fn add(self, rhs: Self) -> Self {
        let mut sum = Polynomial::new();

        for (idx,(l_c, r_c)) in self.iter().zip(rhs.iter()).enumerate() {
            sum[idx]=(l_c+r_c).rem_euclid(CHARACTERISTIC);
        }
        sum
    }
}

impl std::ops::Mul for Polynomial {
    type Output = Polynomial;

    fn mul(self, rhs: Self) -> Self {
        let mut prod: Polynomial = Polynomial::new();

        for (o_idx, l_c) in self.iter().enumerate() {
            for (i_idx, r_c) in rhs.iter().enumerate() {
                if o_idx+i_idx > POL_SIZE && r_c*l_c > 0 {
                    panic!();
                }
                if let Some(mutref) = prod.get_mut(o_idx+i_idx) {
                    *mutref += l_c*r_c;
                };
            }
        }

        for coefficient in prod.iter_mut() {
            *coefficient=coefficient.rem_euclid(CHARACTERISTIC);
        }
        prod
    }
}
