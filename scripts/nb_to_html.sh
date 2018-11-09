#!/usr/bin/env bash
file=$1
filepath="${PWD}/${file}"
filename=$(basename "${file}")
filename="${filename%.*}"

mkdir -p tmp
cd tmp

echo """{%- extends 'full.tpl' -%}
{% block in_prompt %}
{% endblock in_prompt %}
{% block input %}
    {% if cell.metadata.suppress_input is defined and
       cell.metadata.suppress_input == 1 %}
    {% else %}
        {{ super() }}
    {% endif %}
{% endblock input %}
{% block output %}
    {% if cell.metadata.suppress_output is defined and
       cell.metadata.suppress_output == 1 %}
    {% else %}
        {{ super() }}
    {% endif %}
{% endblock output %}

""" > my_html.tpl
ipython nbconvert --quiet --to html \
	--template my_html "${filepath}" $> /dev/null && \
google-chrome "${filename}.html" $> /dev/null &
