# Crack Propagation due to Spatially Varying Temperature Field

## Motivation

Problems related to growth and fracture have been a topic of interest for some time now. Not only these systems are rich in physics but we also get to observe some very beautiful patterns and colors, for example check out this work by [Atis et al](https://journals.aps.org/prx/abstract/10.1103/PhysRevX.9.021058) about how a yeast colony tears itself apart in search for nutrients:

![Atis et al](figures/yeast.gif)

Here's another example from [Plummer et al](https://pubs.rsc.org/en/content/articlelanding/2024/sm/d3sm01470c) in which a hydrogel is allowed to swell in a constrained environment.

![Plummer et al](figures/hydrogelswelling.gif)

Even though the field of fracture mechanics is at least a century old, but still it is very difficult to be confident about the numerical simulations one is running related to it. And the literature related to growth and fracture is not very vast but thankfully, growth strains and thermal strains are implemented on a system in exactly similar fashion (Please note that the discussion here is limited to implementing growth or thermal strains, but not the physics behind how those strains are being produced in the first place. For more information, please check out the chapters 13 and 14 of this book on growth by [Alain Goriely](https://link.springer.com/book/10.1007/978-0-387-87710-5)). 

The deformation gradient $F$ can be decomposed into the product of two tensors, $A$ due to residual stresses and $G$ (or $T$ for thermal strain) due to growth

$$F = AG$$

Where $F = I + \mathr{grad}(u)$ as usual.
