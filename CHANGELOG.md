# Changelog

## [0.4.0](https://github.com/hpc-maths/ponio/compare/v0.3.0...v0.4.0) (2025-10-24)


### Features

* adaptive time step for Strang splitting method ([#53](https://github.com/hpc-maths/ponio/issues/53)) ([07e9c98](https://github.com/hpc-maths/ponio/commit/07e9c9892c750f51f448945ff9bf343aae75a56c))
* add adaptive time step for samurai and ROCK ([#54](https://github.com/hpc-maths/ponio/issues/54)) ([0eb5e33](https://github.com/hpc-maths/ponio/commit/0eb5e33fec845f892877c18c7f1e05ae5d60a9d6))
* add ESDIRK (5,4) a from Kvaerno 2004 ([#89](https://github.com/hpc-maths/ponio/issues/89)) ([9b88dff](https://github.com/hpc-maths/ponio/commit/9b88dff270748cf81dd9198e0998a731f6023b10))
* add informations on time interation in time_iterator ([#59](https://github.com/hpc-maths/ponio/issues/59)) ([6a237f3](https://github.com/hpc-maths/ponio/commit/6a237f3ed53711f15ee0ef33001f138c92098e4b))
* add interface for user to provide an algorithm ([#66](https://github.com/hpc-maths/ponio/issues/66)) ([ea4bc79](https://github.com/hpc-maths/ponio/commit/ea4bc79895f514baeab7b2e6e886af36b98c6f97))
* add option to build ponio offline ([#79](https://github.com/hpc-maths/ponio/issues/79)) ([5bf910f](https://github.com/hpc-maths/ponio/commit/5bf910f7fa94924fdefe79e0740d723f14cee8dc))
* add PIROCK with 3 operators ([#65](https://github.com/hpc-maths/ponio/issues/65)) ([4a2e4a6](https://github.com/hpc-maths/ponio/commit/4a2e4a6eb04721a6af5465e92023b932a45e2ff3))
* add RKL method ([#72](https://github.com/hpc-maths/ponio/issues/72)) ([9172f47](https://github.com/hpc-maths/ponio/commit/9172f47879b580108b1d03861a20c7f1b38e62c7))
* add some SDIRK methods ([#71](https://github.com/hpc-maths/ponio/issues/71)) ([5c05dfa](https://github.com/hpc-maths/ponio/commit/5c05dfafae68d27c74f48e0c7b137ef9e44bfa8d))
* Add two Butcher tableaus ([#74](https://github.com/hpc-maths/ponio/issues/74)) ([3d1804a](https://github.com/hpc-maths/ponio/commit/3d1804a63046a74de4a1242d760f004496504da2))
* change default build to -DBUILD_OFFLINE ([#80](https://github.com/hpc-maths/ponio/issues/80)) ([c1a5744](https://github.com/hpc-maths/ponio/commit/c1a57444a2fed7920ae3077dc3de89a33e724cc9))
* change samurai version in examples ([#86](https://github.com/hpc-maths/ponio/issues/86)) ([10ed306](https://github.com/hpc-maths/ponio/commit/10ed306e9467456248e14ae04e36f4ed59f3ff8c))
* improve adaptive Strang method with estimate of Lipschitz constant ([#83](https://github.com/hpc-maths/ponio/issues/83)) ([b9e56fa](https://github.com/hpc-maths/ponio/commit/b9e56fa5b8f3c6684ad9cf05f23e2d57cd7dd688))
* improve CI, doc and tests ([#55](https://github.com/hpc-maths/ponio/issues/55)) ([79078e5](https://github.com/hpc-maths/ponio/commit/79078e5000604bfb5cea39ea38da35d839a0103d))
* improve PIROCK adaptive time step method ([#76](https://github.com/hpc-maths/ponio/issues/76)) ([19760e7](https://github.com/hpc-maths/ponio/commit/19760e708f17cc9c98212a341f022d30be17b786))


### Bug Fixes

* change iteration_info interface on algorithms ([#70](https://github.com/hpc-maths/ponio/issues/70)) ([f988bc7](https://github.com/hpc-maths/ponio/commit/f988bc738dcf9895da07a944f17a92bd946fec13))
* change version of release-please and add token ([#64](https://github.com/hpc-maths/ponio/issues/64)) ([d183ae5](https://github.com/hpc-maths/ponio/commit/d183ae570f5779781792829e68a7a79dd733fafe))
* embedded ROCK4 ([#69](https://github.com/hpc-maths/ponio/issues/69)) ([7982712](https://github.com/hpc-maths/ponio/commit/7982712156858d96e5a3a9a237310dec211a1f6a))
* fix DIRK methods with samurai ([#88](https://github.com/hpc-maths/ponio/issues/88)) ([be7bce9](https://github.com/hpc-maths/ponio/commit/be7bce9f3a065154c945a1b95744704a69652307))
* fix midnight CI ([#56](https://github.com/hpc-maths/ponio/issues/56)) ([783da0a](https://github.com/hpc-maths/ponio/commit/783da0a6a64e247f0c21097751d956861b9b8291))
* fix pixi version for pages ([#87](https://github.com/hpc-maths/ponio/issues/87)) ([df2c48d](https://github.com/hpc-maths/ponio/commit/df2c48d0dc0b596e5a71974117a86350925b961c))
* remove brackets and fix coefficients ([#85](https://github.com/hpc-maths/ponio/issues/85)) ([dbff331](https://github.com/hpc-maths/ponio/commit/dbff331b97cd8fbf01cdbcd4f6bae764374fa709))
* remove const for samurai solver ([#82](https://github.com/hpc-maths/ponio/issues/82)) ([99a1398](https://github.com/hpc-maths/ponio/commit/99a1398938147efd2f212ac4a659869c331326d9))
* remove first iteration with time_step equals to zero ([#57](https://github.com/hpc-maths/ponio/issues/57)) ([15d096c](https://github.com/hpc-maths/ponio/commit/15d096ca34016c4be401cd2c1623b51c10e71eb6))
* remove warning and error of compilation with NonLinearLocalSolvers ([#61](https://github.com/hpc-maths/ponio/issues/61)) ([f438db5](https://github.com/hpc-maths/ponio/commit/f438db5223dc9f7a6e23b15a0f597adecd1bd0dd))
* restore version because unrelease v0.3.0 ([#63](https://github.com/hpc-maths/ponio/issues/63)) ([17a1b6f](https://github.com/hpc-maths/ponio/commit/17a1b6fda7b96c5b651703784747ddfe76403ca4))
* Shampine's trick with samurai ([#52](https://github.com/hpc-maths/ponio/issues/52)) ([5981006](https://github.com/hpc-maths/ponio/commit/598100653983dee98b4ca994ad6ea887f695c336))


### Performance Improvements

* change given function from f(t, y) to f(t, y, dy) ([#84](https://github.com/hpc-maths/ponio/issues/84)) ([8aeaae1](https://github.com/hpc-maths/ponio/commit/8aeaae195dd16310b99a40a526726701111d53dc))

## [0.3.0](https://github.com/hpc-maths/ponio/compare/v0.2.0...v0.3.0) (2024-05-30)


### Features

* add PIROCK method ([#44](https://github.com/hpc-maths/ponio/issues/44)) ([9df536f](https://github.com/hpc-maths/ponio/commit/9df536f0c17719af8f0a5086ea7e31b293b76111))

## [0.2.0](https://github.com/hpc-maths/ponio/compare/v0.1.0...v0.2.0) (2024-03-18)


### Features

* rename solver into ponio ([#47](https://github.com/hpc-maths/ponio/issues/47)) ([b9419c8](https://github.com/hpc-maths/ponio/commit/b9419c8a0e907f97e54c24da53b7d08d3bd960cf))


### Bug Fixes

* cmake ([#46](https://github.com/hpc-maths/ponio/issues/46)) ([810e93e](https://github.com/hpc-maths/ponio/commit/810e93e767cb358e5e1fa337e654d5e9b64397a4))
