// Spoken narration for the two films, one string per caption beat of each chapter (keyed by chapter
// id), plus the title cards. The text is written to be read aloud by a speech synthesizer, so all
// symbols and numbers are spelled out; the on-screen captions in chapters.js stay authoritative.
(function (root) {
  'use strict';
  root.STBV.NARRATION = {
    _film: 'The spacetime bootstrap. How the spin-one Haldane gap is proved, starting from imaginary time.',
    _open: 'Part two of the spacetime bootstrap: open chains, a chain with two ends.',

    gap: [
      'A chain of L spins has a ladder of energy levels. The gap, Delta, is the first step up the ladder.',
      'In nineteen eighty-three, Haldane predicted that for spin one the gap stays open as the chain grows, while for spin one-half it closes.',
      'Our goal is one number c, greater than zero, such that the gap Delta of L is at least c, for every ring length L.'],
    rate: [
      'At inverse temperature beta, level m has probability e to the minus beta E sub m, divided by Z.',
      'Cooling down empties the excited levels. The ratio of the first two weights is exactly e to the minus beta Delta.',
      'So it is enough to show that the excited weight is at most e to the minus c beta, on every ring. Then the gap is at least c.'],
    imag: [
      'A quantum state is a list of numbers, c sub m, one for each energy level. Draw each one as an arrow: its length is the size of c sub m, and its angle is the phase. The bars below are the probabilities, the squared lengths.',
      'In real time, e to the minus i H t, arrow m turns at speed E sub m. The lengths never change, so the probabilities never change. Nothing settles.',
      'Now replace t by minus i tau. The same formula gives e to the minus tau H, with no turning at all. Arrow m shrinks by the factor e to the minus tau E sub m, and high energies shrink fastest.',
      'Rescale so that the probabilities add up to one. The ground state takes over, and the first excited level loses out like e to the minus tau Delta. The gap is the speed.',
      'At tau equal to beta, the shrink factors are exactly the Boltzmann weights, e to the minus beta times E. Low temperature and long imaginary time are the same thing.'],
    paths: [
      'Start with how matrices multiply. The entry of a product, row i, column k, is a sum over a middle index j: to go from i to k, pass through some middle state j, multiply the two entries you use on the way, and add over every j.',
      'With numbers. Take the matrix two, one, one, three, and square it. The top left entry of the square has two routes: through state one, two times two, which is four, and through state two, one times one, which is one. Four plus one is five.',
      'Now repeat for M factors. An entry of a matrix raised to the power M is a sum over paths through M steps. Each path picks one entry of the matrix at each step, and its weight is the product of those entries. That is why weights multiply along a path.',
      'Take one spin that can flip, with H equal to minus sigma x. For a tiny step d tau, the step matrix, e to the minus d tau H, is about one on the diagonal and d tau off it. Staying costs one, and flipping costs d tau. Then e to the minus beta H is the step matrix to the power M, with M equal to beta over d tau.',
      'Make the steps tiny. A path becomes a world line that mostly stays put, and flips a few times. With k flips its weight is d tau to the k, and there are about beta over d tau, to the k, over k factorial, places to put them. In total, beta to the k over k factorial.',
      'The trace keeps only the paths that end where they started. So k is even, and the time axis closes into a circle.',
      'Add up all the paths: two beta to the k over k factorial, summed over even k, gives two cosh beta, the exact answer. This is the Euclidean path integral: quantum mechanics as a sum over histories.'],
    torus: [
      'The partition function Z is the trace of e to the minus beta H. Draw the chain across, and imaginary time going up. That makes a sheet.',
      'The trace makes the last time slice equal to the first, so imaginary time is a circle.',
      'The chain is a ring, so space is a circle too. Every quantity in the proof lives on this torus.'],
    ising: [
      'Before the quantum chain, a classical warm-up: the Eesing ring. Each site carries a spin, plus or minus one, and each pair of neighbors contributes a weight: e to the K if they agree, and e to the minus K if they disagree.',
      'The weight of a whole configuration is the product of its bond weights. Here that is e to the minus two K. The partition function Z adds up the weights of all two to the L configurations.',
      'Store the bond weight in a two by two matrix T, whose entry for spins s and s prime is e to the K s s prime. Draw the matrix as a box with two legs, one for each index. The spins now sit on the legs between boxes.',
      'A leg joining two boxes means: sum over that index. Summing the spin s two, shared by the first two boxes, gives exactly the matrix product: the entry of T squared from s one to s three. This is what it means for T to pass the spin on to the next site.',
      'Now join every leg, one after another. The row of boxes becomes T to the L. The last leg wraps around the ring back to the first site, and joining those two ends is a trace. So Z is the trace of T to the L.'],
    chain: [
      'Now the quantum chain. As for the single spin, cut beta into tiny steps, d tau. Then e to the minus beta H is about one minus d tau H, to the power M. In each step, either nothing happens, with weight one, or one term of minus H acts, with weight d tau.',
      'For the Heisenberg chain, minus S j dot S j plus one is a sum of three terms, G alpha on one site times G alpha on the next, with alpha equal to x, y, or z. So one term acting means: pick a bond and a type, x, y, or z, and apply G of that type to both sites of the bond.',
      'What are these G matrices? In the basis x, y, z, they are three real antisymmetric three by three matrices. G z turns x into minus y, turns y into x, and kills z: a quarter turn about the z axis. For example, G z times G z is diagonal minus one, minus one, zero, whose trace is minus two.',
      'Each such term is an event: a rung across one bond, at one imaginary time, labeled by its type. A history is the set of all events. Each event carries weight d tau, exactly like a flip of the single spin.',
      'Follow site three. Its events arrive through its left bond, the gold list, and through its right bond, the blue list: z at time zero point two, y at zero point four, and x at zero point six. Reading upward, its matrices multiply in time order, with later ones on the left.',
      'As a tensor, call this W. It has two time legs, s prime at tau equal to zero and s at tau equal to beta, and W is the matrix element of G x G y G z between them. Imaginary time is a circle, so join the two time legs. That sets s equal to s prime, a trace. So X is the trace of G x G y G z, which is one.',
      'Other cases: no events gives the trace of the identity, three. One z event on each side gives the trace of G z times G z, minus two. The weight of a whole history is a trace over the whole chain, and it factorizes into one such trace per site.',
      'So a history weighs a product of X\'s, one per site, each taking the histories of the two bonds next to it. In pictures: boxes X on the sites, joined by legs that carry bond histories. It is the Eesing diagram with T replaced by X. So Z is the trace of X to the L. This X is the transfer matrix of the proof.'],
    slice: [
      'Slice across: each row is a step of time evolution. Summing over energies gives Z as a sum of e to the minus beta times E.',
      'Slice up instead: each column is one site, carried along the ring by a transfer matrix X. Then Z is the trace of X to the L, a sum of lambda to the L.',
      'The same number, as two sums. On a long ring the largest lambda wins, just as at low temperature the lowest energy wins.'],
    purity: [
      'A state with probabilities p one, p two, and so on, has purity P, the sum of p squared. Draw each p as a square: the purity is their total area.',
      'Every square fits under the height of the largest one, p zero. So the purity is at most p zero.',
      'A purity near one therefore forces a single state to carry almost all of the probability.'],
    two: [
      'Draw the torus as a network: one transfer matrix X per site, joined leg to leg around the ring. Joining all the legs gives Z n, the trace of X to the n.',
      'Now cut the torus along imaginary time, between two sites. That cuts one leg of the ring. What is left is the matrix X to the n, with two open legs. Divide by Z n so that its trace is one. That is sigma, with eigenvalues p sub i equal to lambda sub i to the n over Z n.',
      'The purity of sigma joins two copies leg to leg and closes the ring. The trace of sigma squared is the trace of X to the two n over Z n squared, which is Z two n over Z n squared: a ring twice as long. This is the spatial purity.',
      'Or cut along space instead, at one moment of imaginary time, through every site. What is left is e to the minus beta H, with open legs above and below. Divide by the partition function and you get the thermal state rho, whose eigenvalues are the Boltzmann probabilities.',
      'The purity of rho stacks two copies and closes time. The trace of rho squared is Z n at two beta, over Z n at beta squared: a torus twice as tall. This is the thermal purity.',
      'Both purities are ratios of partition functions of small tori, numbers that a computer can evaluate.'],
    square: [
      'Square every probability, and renormalize. The leader pulls ahead.',
      'Why is this gluing? Sigma squared over its trace is X to the two n over Z two n, the state of the ring twice as long. Its probabilities are p sub i squared, over the sum of the squares. So gluing two copies squares every probability, and the same holds in time.',
      'A short lemma bounds what squaring does to the defect epsilon, one minus the purity. The new defect is at most epsilon squared over two times one minus epsilon squared. That is about epsilon squared over two.',
      'Start below one half, and repeated squaring drives the defect to zero: the cobweb drops off the diagonal.'],
    couple: [
      'Squaring grows the torus in one direction at a time. We need both.',
      'First, bookkeeping: every gluing multiplies by a purity. Side by side, Z two n is the spatial purity times Z n squared. Stacked, Z n at two beta is the thermal purity times Z n squared. These are just the two definitions, rearranged.',
      'Start from four copies of the torus, and build one torus twice as long and twice as tall. Rows first: glue side by side, which multiplies by the spatial purity squared, then stack, which multiplies by the thermal purity at length two n.',
      'Or columns first: stack, which multiplies by the thermal purity squared, then glue side by side, which multiplies by the spatial purity at two beta. Both routes end at the same torus.',
      'So the two products are equal. The spatial purity at two beta, times the thermal purity squared, equals the spatial purity squared, times the thermal purity at length two n. This coupling identity ties the spatial and thermal purities together.'],
    rpwhy: [
      'The squaring step reads p sub i, equal to lambda sub i to the n, over Z n, as probabilities. Let us test that on a two by two transfer matrix, with actual numbers.',
      'First, a symmetric matrix. Its eigenvalues, one point zero six and zero point three four, are real. Compute Z n, the trace of X to the n, and the spatial purity, Z two n over Z n squared, for n from one to eight.',
      'Every purity lies between zero and one, and climbs to one. This is the behavior the bootstrap relies on.',
      'Now a matrix that is not symmetric: a small rotation. Its eigenvalues, zero point nine plus or minus zero point four i, are complex.',
      'Now Z n changes sign, and the so-called purity becomes negative, or huge: minus forty-seven at n equal to four. These numbers are no longer probabilities, so squaring proves nothing. The proof needs X to be symmetric.'],
    rpwhere: [
      'In the proof, X is a matrix. Its row index is its left leg, the history of the bond on the site\'s left, and its column index is its right leg, the history of the bond on its right. Symmetric means that swapping the two gives the same number: mirroring the box, left to right, changes nothing.',
      'Let us test the mirror on one site. An event\'s side is the bond it arrives through: the left bond, shared with the left neighbor, or the right bond, shared with the right neighbor. For example, G z at time zero point two and G x at zero point six arrive on the left, and G y at zero point four on the right.',
      'For the Heisenberg chain, an event puts the same matrix on both sites of its bond, so the site gets the same G whichever side the event comes from. In time order it sees G z, then G y, then G x. So X is the trace of G x G y G z, which is one.',
      'Now mirror it: the left events move to the right, and the right event to the left. The site sees the same three matrices at the same times, so X is one again. The Heisenberg X is symmetric.',
      'Now a counter-example: a bond whose event puts G z on its left site, but P, the projector onto x, on its right site. A site is the right end of its left bond, so it receives P from there. And it is the left end of its right bond, so it receives G z from there.',
      'Two events through the left bond give the trace of P times P, which is one. The same two events through the right bond give the trace of G z times G z, minus two. So X is not symmetric. Reflection positivity is the condition that rules this out: the bond must look the same from both sides, with positive weights.'],
    rptest: [
      'Here is a test that any bond can be put through. For the Heisenberg chain the bond is h, equal to S j dot S j plus one, so minus h is the sum of G alpha tensor G alpha. As a tensor, minus h has four legs: i and k on the left site, its output and input, and j and l on the right site. Each entry is a sum over alpha of the two G entries.',
      'The usual way to read it as a matrix: the rows are the outputs, i and j, and the columns are the inputs, k and l. That is this nine by nine matrix, minus h for two spin-one sites. For example, the entry from y y to x x is G z x y times G z x y, which is one.',
      'Realignment means reading the same eighty-one numbers with the legs grouped differently: the rows are the left site\'s legs, i and k, and the columns are the right site\'s legs, j and l. Nothing is computed. Each entry just moves to its new place. The entry from y y to x x goes to row x y, column x y.',
      'Why group the legs this way? Minus h is two boxes, G alpha on each site, joined by a summed leg alpha. Split between the sites, each half is a vector, v alpha, made from the entries of G alpha. The joined leg is a matrix product, so the realigned minus h is the sum of v alpha times v alpha transpose: a sum of squares.',
      'A sum of squares has no negative eigenvalues. Here they are two, two, two, and then zeros. The antiferromagnet passes.',
      'For the ferromagnet, h is minus S dot S, so minus h is minus the sum of G alpha tensor G alpha, and the realigned matrix has eigenvalues minus two. That is not a sum of squares, and not reflection positive. In general, if the realigned minus h has no negative eigenvalues, up to an on-site shift, the bond is reflection positive.'],
    double: [
      'Start with both purities at least one minus epsilon, at size n and inverse temperature beta: the spatial purity at length n, and the thermal purity at length two n.',
      'Move one: couple. The coupling identity gives the spatial purity at two beta as the spatial purity squared, times the thermal purity at two n, divided by the thermal purity at n, squared. Purities are at most one, so dividing only helps. The result is at least one minus epsilon, cubed: a defect of about three epsilon.',
      'Move two: square in space. Gluing two rings sends the defect three epsilon to about three epsilon squared over two, which is four and a half epsilon squared, below five epsilon squared. That bounds the spatial purity at two n and two beta.',
      'Move three: couple again, now at size two n. The same identity, solved for the thermal purity at four n, gives at least one minus five epsilon squared, times one minus epsilon, squared: a defect of about two epsilon.',
      'Move four: square in time. Stacking sends the defect two epsilon to about two epsilon squared, again below five epsilon squared. That bounds the thermal purity at four n and two beta.',
      'Out comes the same pair of bounds at two n and two beta, with epsilon replaced by five epsilon squared. The step repeats.'],
    iterate: [
      'Each step doubles beta, and roughly squares epsilon. Plot the log of epsilon against beta.',
      'The points fall on a straight line: epsilon decays exponentially in beta.',
      'Why is that a gap? A thermal purity of at least one minus epsilon forces the ground state to hold probability at least one minus epsilon, because the purity is at most the largest probability. Then the first excited level has probability at most epsilon, and the ratio of the two is e to the minus beta Delta.',
      'So e to the minus beta Delta is at most epsilon over one minus epsilon, and Delta is at least the log of one minus epsilon over epsilon, divided by beta: the slope of the line. For spin one, the gap is greater than log twenty over seven hundred eighty-four, about zero point zero zero three eight, for every even ring length of at least two thousand three hundred four.'],
    seedtoy: [
      'The doubling step needs a starting point: a purity close to one at some size n zero. The spatial purity is the sum of p sub i squared, where p sub i is lambda sub i to the n, over Z n. So it depends on the eigenvalues, lambda, of the transfer matrix X.',
      'But X acts on whole histories, so it is far too big to diagonalize. What a computer can get exactly is Z n, the trace of X to the n, which is the sum of lambda to the n. That is just the partition function of a small ring of n sites. The question is: can a few such numbers pin down the purity?',
      'Let us try it on a small made-up example. Pretend that X has four eigenvalues: one, zero point six, minus zero point five, and zero point four. And play by the rules: we may use only the traces Z n, never the eigenvalues themselves.',
      'Step one: moments. Z six, the partition function of the six-site ring, is the sum of lambda to the sixth. The ring is small enough to compute it exactly: one point zero six six four. Every term is at least zero, so no single eigenvalue can exceed one point zero six six four to the one sixth, which is one point zero one one, in size.',
      'Step two: a trial state. Z twenty-four, the sum of lambda to the twenty-fourth, belongs to the twenty-four-site ring, which is too big to diagonalize. But Z is a sum of positive Boltzmann factors, so it is at least the ground-state term, e to the minus beta E zero. And for any trial state, the ground energy is at most the trial state\'s energy, which is cheap to compute for a matrix-product state. If that energy is below zero, after shifting all energies by a constant per site, which cancels in the purities, then the ground-state factor exceeds one, and Z twenty-four is at least one.',
      'Now suppose every eigenvalue were at most zero point nine seven in size. Then each lambda to the sixth is at most zero point eight three three, and together they add up to one point zero six six four. Lambda to the twenty-fourth is lambda to the sixth, to the fourth, and fourth powers reward concentration. So the worst case puts one of them right at the cap, and the remaining zero point two three three in one more. Even then, the twenty-fourth moment is at most zero point eight three three to the fourth, plus zero point two three three to the fourth: zero point four eight four, less than one. A contradiction.',
      'So one eigenvalue, lambda zero, exceeds zero point nine seven, and its sixth power is at least zero point eight three three. The others share what is left of the sixth moment: at most one point zero six six four minus zero point eight three three, which is zero point two three three. For numbers that are at least zero, a sum of fourth powers is at most the fourth power of the sum. So the others carry at most zero point two three three to the fourth, or zero point zero zero three, of the twenty-fourth moment.',
      'Step three: the purity. Recall that the spatial purity at twenty-four is Z forty-eight over Z twenty-four squared: the sum of lambda to the forty-eighth, over the square of the sum of lambda to the twenty-fourth. On top, every term is at least zero, so keeping only lambda zero makes it smaller. Below, the sum is lambda zero to the twenty-fourth plus the rest, and the rest is at most zero point zero zero three.',
      'Putting these together, the purity is at least one over the square of one plus zero point zero zero three over lambda zero to the twenty-fourth. Since lambda zero exceeds zero point nine seven, that is at least zero point nine eight eight. A certified starting purity, with defect about zero point zero one two, found without ever knowing the spectrum. The true value is zero point nine nine nine nine nine.',
      'The real certificate follows these same three steps on the spin-one chain, with exact integers throughout. The next chapter draws its ingredients.'],
    seed: [
      'The real certificate, built from the inputs of OpenAI\'s twenty twenty-six proof, runs the same three steps on the spin-one chain. Step one needs exact moments, Z n, the trace of X to the n: the closed ring of n boxes X. That is the trace of e to the minus beta H on a ring of up to twelve sites, evaluated in exact integer arithmetic.',
      'A sharper version inserts a symmetry operation, g, into the ring before closing it. These twisted traces separate the eigenvalues by symmetry sector, so each sector gets its own moment bound.',
      'Step two needs a trial state on seventy-two sites: a matrix-product state. It has one box, A, per site, with a physical leg for the spin, and two virtual legs of dimension twenty-six, joined around the ring.',
      'Its norm glues the state to its mirror image along the physical legs. That makes a double-layer ring: the trace of T to the seventy-second power, where T joins one column of the two layers. Its energy is seventy-two times the same ring with one h inserted. Both are rings of small matrices, and they are contracted exactly.',
      'Why does that give Z seventy-two greater than one? The energies carry a constant shift, kappa per site, which cancels in every purity. So Z seventy-two is a sum of positive Boltzmann factors, and it is at least the ground-state factor. The variational principle puts the ground energy below the trial energy, and the trial energy is below minus seventy-two kappa. So that single factor already exceeds one.',
      'Then, as in the made-up example: Z seventy-two is the sum of lambda to the seventy-second, it exceeds one, and the moment bounds cap every eigenvalue. Together they force one eigenvalue above zero point nine nine nine. That certifies both purities at ring length thirty-six.',
      'The doubling step, run six times, carries them to length two thousand three hundred four, at beta seven hundred eighty-four, with both defects below zero point zero one.'],

    owhy: [
      'The ring proof needs a ring. Z, the trace of X to the n, joins the last leg back to the first. But the boundary-field theorem is about a chain with two ends. Cut the ring open, and there is no loop left to close.',
      'A picture of spin one helps here: the valence-bond solid, or V B S, of Affleck, Kennedy, Lieb and Tasaki. Write each spin one as two spin one-halves, the two dots in one oval, combined symmetrically.',
      'Neighbors pair one half each into a singlet, with total spin zero, shown in blue. On a ring, every half finds a partner, and nothing is left free.',
      'Now cut the ring open. The two outermost halves lose their partners. A free spin one-half sits at each end, which gives two times two, four, low states: a singlet and a triplet, whose splitting shrinks as the chain grows.'],
    oedge: [
      'Here are real numbers with no field, from exact diagonalization, for n from four to twelve. The end states sit low. For even n, a triplet sits just above the ground state, at zero point five one, zero point three one, zero point two, zero point one four, and zero point one. For odd n, the ground state itself is a triplet. There is no uniform gap.',
      'An end field, h equal to three fifths along z, pins the end spins and splits the four states. For odd n the ground state is alone, with gaps zero point seven two, zero point six four, zero point five eight, and zero point five five, leveling off. Even n are worse: the left and right ends prefer opposite directions.',
      'That is the chain of the theorem: H is the sum of S dot S, minus h times the z spins at the two ends, with h equal to three fifths, and n odd.'],
    oising: [
      'To see how an end is closed, go back to the classical Eesing chain, now open. Each of its n minus one bonds still gives one matrix K, joined leg to leg. But the first and last spins have an outer leg with nothing attached.',
      'An end field gives each end spin an extra weight, e to the beta h s: a vector u, with one entry for each spin value. Attach it to the free outer leg as a cap.',
      'Summing every spin now gives Z as u, then K, K, and so on, then u: u transpose, K to the n minus one, u. A ring closes with a trace; an open chain closes with two vectors.',
      'In terms of eigenvalues, u transpose K to the m, u, is the sum of lambda i to the m, times the squared overlap of u with the eigenvector phi i. A ring weighs every eigenvalue by one. An open chain weighs it by that overlap, and misses lambda zero entirely if u is orthogonal to its eigenvector.'],
    ocyl: [
      'Now the quantum chain. For a ring, space and imaginary time are both circles: the torus.',
      'Open the chain. Imaginary time is still a circle, because Z is still a trace over spin states, but space now has two ends. The spacetime is a cylinder.',
      'Inside, nothing changes: the same box X sits at every site. Only the two ends need closing, as the vectors u did for Eesing. Those caps, v bar and v, are the boundary vectors.'],
    ovec: [
      'What closes an end? The end site has a bond on its inner side only. On its outer side acts the end field, minus h S z. It enters e to the minus beta H just like a bond term, so it fires field events at random imaginary times.',
      'A bond event puts G alpha on both of its sites. A field event is half a bond: it puts G z on the end site alone, with a factor i h, since S z is minus i G z. Here are histories with zero, one, two, and three field events.',
      'The field events form a history, just like the missing bond\'s. So the end box\'s outer leg is capped by a vector of histories: v of a history is i h to the power k when all k events are of type z, and zero otherwise. Then the open partition function is v bar, X to the n, v, just like u transpose K to the m, u.',
      'The left cap uses the complex conjugate, v bar, whose events carry minus i h: a field of the opposite sign. So v bar, X to the n, v is the chain with opposite end fields. Its two caps are the same vector, so its even powers are squared norms, like the ring\'s.',
      'The equal-field chain of the theorem is reached by U, the rotation by pi about the x axis, which turns S z into minus S z. It flips the sign of one end field, giving v bar, X to the n, U, v.'],
    over: [
      'Normalize so that the top eigenvalue is one: T is X over lambda zero. The ring estimates show that lambda zero is simple and every other eigenvalue is smaller. So T is the projector onto psi, its top eigenvector, plus a remainder R, with norm less than one.',
      'Apply T to the boundary vector, again and again. The part along psi stays; the rest shrinks by the norm of R each time.',
      'So the open partition function, divided by lambda zero to the m, is B m, equal to c plus D m: a fixed overlap c, the squared overlap of psi with v, as for Eesing, plus a remainder D m that dies away with the length m.',
      'The boundary defect, e, equal to D m over B m, is the share not carried by psi. Any e less than one proves that c is positive: a long chain does see the ground state.'],
    ofold: [
      'The remainder D m is v, R to the m, v, for even m. Write R to the m as R to the m over two, times R to the m over two, and fold the strip at its middle.',
      'R is symmetric, because X is, by reflection positivity. So the two halves are mirror images, and D m is the squared norm of R to the m over two, applied to v. It is never negative.',
      'Now move k over two factors of R onto each end, with w equal to R to the k over two, applied to v, so that D k is the squared norm of w. The middle, R to the m minus k, has norm at most the norm of R to the m minus k. So D m is at most that, times D k.',
      'Longer chains shrink the remainder geometrically, while the overlap c stays put. That is what lets the boundary defect improve with length.'],
    omoves: [
      'Here is the plan. For a long open chain, the partition function is about c times lambda zero to the n: one bulk factor per site, and one boundary factor, c. So its thermal purity is about a bulk ratio of eigenvalues, to the n, times a boundary ratio of overlaps.',
      'The bulk factor is a ring quantity, already controlled by the ring proof. The boundary factor does not depend on n, so it cancels between lengths M and two M. That works only while the leading term dominates, so the boundary defect is tracked too.',
      'At each scale j, track two numbers: the boundary defect e, and the open-chain thermal purity at length M j. Both defects start at most x.',
      'Move one: the overlap recurrence. The remainder suppression of a longer chain, at the next temperature, squeezes the boundary defect to at most ten x squared.',
      'Move two: double the length. The purity at two M j is at least the product of one minus the old open defect, one minus the ring defects, and one minus the new boundary defect: about one minus two x.',
      'Move three: square the Gibbs state, which doubles beta and needs no ring. The defect drops to ten x squared.',
      'Out comes the same pair at the next scale, with x replaced by ten x squared. The step repeats, just as on the ring.'],
    onum: [
      'Here is the note\'s first step, from scale two to scale three, with every number. Given: a boundary defect of zero point one four, an open-chain thermal defect, eta, of zero point two one, and periodic defects epsilon two, zero point zero five two, and epsilon three, six epsilon two squared, zero point zero one six two.',
      'Move one. The suppression factor q is zero point one four two. The overlap recurrence bounds the new boundary defect by q over one minus q, times a bracket built from the old defects. That is zero point one six five five times zero point seven one one five, or zero point one one eight, below the table entry of zero point one two.',
      'Moves two and three. With the boundary defect rounded up to zero point one two, the purity after doubling is at least zero point seven nine, times zero point nine four eight, times zero point nine eight four, times zero point eight eight, which is zero point six four eight. Squaring turns the defect, zero point three five two, into zero point one four seven, below zero point one five five.',
      'Both defects went down. Two more such steps reach zero point zero zero one. After that, the clean rule, x goes to ten x squared, takes over.'],
    ostart: [
      'The open-chain iteration starts at beta equal to forty-nine, from four certified numbers: two periodic purities at length one hundred twenty, and two open-chain bounds from the boundary-field paper.',
      'Three rounds of the moves take the boundary defect from zero point one four, and the thermal defect from zero point two one, down to zero point zero zero one, at beta three hundred ninety-two and length one thousand nine hundred twenty.',
      'From there, x goes to ten x squared: ten to the minus two, then minus three, minus five, minus nine, minus seventeen. Doubly exponential, exactly as on the ring.'],
    ogap: [
      'The theorem is for odd n. R can have negative eigenvalues, so for odd n the remainder need not be positive, and the fold fails. Instead, put M over two factors of R on each end, so that w is R to the M over two, applied to v.',
      'The middle part, R to the n minus M, times U, has norm at most one. So the remainder is at most the squared norm of w, which is D of M, whatever its sign.',
      'Every n of at least one thousand nine hundred twenty lies in a window from M j to two M j, with beta j equal to three hundred ninety-two times two to the j minus five. In that window, the thermal purity is at least one minus four x j, which is more than one half.',
      'So the ground state is unique, and e to the minus beta j Delta is less than ten to the minus two to the j minus five. Dividing by beta j, the gap is greater than log ten over three hundred ninety-two, about zero point zero zero five nine, for every odd n of at least one thousand nine hundred twenty-one.'],
  };
})(typeof window !== 'undefined' ? window : globalThis);
