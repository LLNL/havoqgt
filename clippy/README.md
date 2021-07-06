# CLIPPy - Command Line Interface Plus Python
 ````
 ╭────────────────────────────────────╮
 │ It looks like you want to use HPC. │ 
 │ Would you like help with that?     │
 ╰────────────────────────────────────╯
  ╲
   ╲
    ╭──╮  
    ⊙ ⊙│╭
    ││ ││
    │╰─╯│
    ╰───╯
````

## Overview

Clippy (CLI + PYthon) is a Python language interface to HPC resources. Precompiled binaries
that execute on HPC systems are exposed as methods to a dynamically-created `Clippy` Python
object, where they present a familiar interface to researchers, data scientists, and others.
Clippy allows these users to interact with HPC resources in an easy, straightforward
environment – at the REPL, for example, or within a notebook – without the need to learn
complex HPC behavior and arcane job submission commands.


## Building C++ Examples & Maptrix on LC

```console
$ . /usr/workspace/llamag/spack/share/spack/setup-env.sh
$ spack load gcc
$ spack load boost
$ cd cpp
$ mkdir build
$ cd build
$ cmake ../
$ make
$ cd ../.. #back to root project directory 
```

## Running Current Examples
```python
$ ipython3-3.8.2

In [1]: from clippy import Clippy

 ╭────────────────────────────────────╮
 │ It looks like you want to use HPC. │ 
 │ Would you like help with that?     │
 ╰────────────────────────────────────╯
  ╲
   ╲
    ╭──╮  
    ⊙ ⊙│╭
    ││ ││
    │╰─╯│
    ╰───╯

In [2]: c = Clippy({'examples': 'clippy/cpp/build/examples'}, cmd_prefix='', loglevel=0)

In [3]: c.examples.howdy('Seth')  # can also use a named arg: c.examples.howdy(name='Seth')
Out[3]: 'Howdy, Seth!'

In [4]: c.examples.sum(1, 2)
Out[4]: 3.0

In [5]: c.examples.sort_edges([(5,5),(3,5),(2,2),(0,0)])
Out[5]: [[0, 0], [2, 2], [3, 5], [5, 5]]

In [6]: c.examples.sort_edges([(5,5),(3,5),(2,2),(0,0)], reverse=True)
Out[6]: [[5, 5], [3, 5], [2, 2], [0, 0]]

In [7]: c.examples.sort_strings(['zulu','yankee','whiskey','uniform','romeo','mike','kilo','foxtrot','delta','alfa'])
Out[7]: 
['alfa',
 'delta',
 'foxtrot',
 'kilo',
 'mike',
 'romeo',
 'uniform',
 'whiskey',
 'yankee',
 'zulu']

In [8]: c.examples.sort_strings(['zulu','yankee','whiskey','uniform','romeo','mike','kilo','foxtrot','delta','alfa'], reverse=True)
Out[8]: 
['zulu',
 'yankee',
 'whiskey',
 'uniform',
 'romeo',
 'mike',
 'kilo',
 'foxtrot',
 'delta',
 'alfa']

In [9]: c.examples.grumpy()
---------------------------------------------------------------------------
ClippyBackendError                        Traceback (most recent call last)
<ipython-input-8-bf61ad375b31> in <module>
----> 1 c.grumpy()

~/clippy/clippy.py in fn(self, *args, **kwargs)
    125             #
    126             # send_dict['args'] = kwargs
--> 127             j = capself.session.exec(capself.name, kwargs)
    128             return j
    129 

~/clippy/clippy.py in exec(self, cmd, submission_dict)
    196         self.logger.debug(f'run(): result = {p}')
    197         if p.returncode != 0:
--> 198             raise ClippyBackendError(p.stderr)
    199 
    200         return json.loads(p.stdout)

ClippyBackendError: terminate called after throwing an instance of 'std::runtime_error'
  what():  I'm Grumpy!
```

## Authors
- Seth Bromberger (seth at llnl dot gov)
- Roger Pearce (rpearce at llnl dot gov)


## License
CLIPPy is distributed under the MIT license.

See [LICENSE-MIT](LICENSE-MIT), [NOTICE](NOTICE), and [COPYRIGHT](COPYRIGHT) for
details.

SPDX-License-Identifier: MIT

## Release
LLNL-CODE-XXXXXX
