#!/bin/bash --login

for i in a b c d
do
	echo $i
	if [[ $i == a ]]
	then
		cp settings.cfg.template settings.cfg.a
		perl -i -pe 's/template_tau/0.0/g' settings.cfg.a
		perl -i -pe 's/template_diffconst/0.536/g' settings.cfg.a
		perl -i -pe 's/template_dirstring/a/g' settings.cfg.a
		./bin/barrio.x settings.cfg.a
		cp ./plt/plot_fields.py.template ./plt/plot_fields.py
		perl -i -pe 's/template_dirstring/a/g' ./plt/plot_fields.py
		./plt/plot_fields.py
	fi
	if [[ $i == b ]]
	then
		cp settings.cfg.template settings.cfg.b
		perl -i -pe 's/template_tau/0.4/g' settings.cfg.b
		perl -i -pe 's/template_diffconst/0.536/g' settings.cfg.b
		perl -i -pe 's/template_dirstring/b/g' settings.cfg.b
		./bin/barrio.x settings.cfg.b
		cp ./plt/plot_fields.py.template ./plt/plot_fields.py
		perl -i -pe 's/template_dirstring/b/g' ./plt/plot_fields.py
		./plt/plot_fields.py
	fi
	if [[ $i == c ]]
	then
		cp settings.cfg.template settings.cfg.c
		perl -i -pe 's/template_tau/0.0/g' settings.cfg.c
		perl -i -pe 's/template_diffconst/0.516/g' settings.cfg.c
		perl -i -pe 's/template_dirstring/c/g' settings.cfg.c
		./bin/barrio.x settings.cfg.c
		cp ./plt/plot_fields.py.template ./plt/plot_fields.py
		perl -i -pe 's/template_dirstring/c/g' ./plt/plot_fields.py
		./plt/plot_fields.py
	fi
	if [[ $i == d ]]
	then
		cp settings.cfg.template settings.cfg.d
		perl -i -pe 's/template_tau/0.4/g' settings.cfg.d
		perl -i -pe 's/template_diffconst/0.516/g' settings.cfg.d
		perl -i -pe 's/template_dirstring/d/g' settings.cfg.d
		./bin/barrio.x settings.cfg.d
		cp ./plt/plot_fields.py.template ./plt/plot_fields.py
		perl -i -pe 's/template_dirstring/d/g' ./plt/plot_fields.py
		./plt/plot_fields.py
	fi
done
