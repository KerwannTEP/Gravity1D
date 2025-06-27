const G = 1.0
const M = 1.0
const L = 1.0
const alpha = 2*L/pi

const kmax = 10
const Jmax = 5.0

const eps_im = 0.01
const nbJ = max(100, floor(Int64, Jmax/(0.2*eps_im)))
const dJ = Jmax/nbJ

const nbp = kmax*nbJ