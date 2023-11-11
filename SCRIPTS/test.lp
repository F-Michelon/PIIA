% Données en entrée
atom(a;b;c;d;e).

% Le script python, qui va définir la fonction gen_id
#script(python)
import itertools
gen = itertools.count(1)

def gen_id():
    "retourne un nouvel identifiant"
    return next(gen)
#end.

% On donne à chaque atome son identifiant, en utilisant la fonction gen_id.
uid(V,@gen_id()):- atom(V).  % le préfix @ permet d'indiquer qu'il s'agit d'une fonction python
#show uid/2.