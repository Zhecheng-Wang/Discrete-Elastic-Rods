# Discrete Elastic Rods Based Simulations
Computational Science and Engineering Certificate Independent Research Project (Fall 2021).

This implementation adapted time-parallel frame in [[Bergou *et al.*, 2010]](http://www.cs.columbia.edu/cg/pdfs/171-threads.pdf) instead of space-parallel frame in [[Bergou *et al.*, 2008]](http://www.cs.columbia.edu/cg/pdfs/143-rods.pdf).

---
### Author
Zhecheng Wang

### Advisor
Prof. Etienne Vouga

---
### Goals
- [X] Implement Discrete Elastic Rods.
- [X] Test the implementation.
- [ ] Application in [insert trade secret here?].

---
### Installation
Clone from this repo

    git clone https://github.com/Zhecheng-Wang/Discrete-Elastic-Rods.git

Initialize build folder and compile the code

    mkdir build
    cd build
    cmake ..
    make

To run the program, run ``discrete_elastic_rods`` in the ``build`` folder.

---
### Executable
#### Arguments
- **[-driver]**
default 0

- **[-test]**
default 0
#### Usage Example
To run test 1 in driver 0, run
``./discrete_elastic_rods -driver 0 -test 1``.
#### Driver Details
- <details>
    <summary> <b>[Driver 0]</b> Discrete Elastic Rods</summary>
    <br> <ul>
          <li><b>[Test 0]</b> testing stretching and bending energy </li>
          <li><b>[Test 1]</b> testing twisting energy </li>
        </ul>
  </details>


---
### References
**Papers:** [[Bergou *et al.*, 2008]](http://www.cs.columbia.edu/cg/pdfs/143-rods.pdf),
[[Bergou *et al.*, 2010]](http://www.cs.columbia.edu/cg/pdfs/171-threads.pdf),
[[Panetta *et al.*, 2019]](https://julianpanetta.com/pdf/xshell.pdf)