from flask import Blueprint, request, render_template, \
                  flash, g, session, redirect, url_for

Vizpymzml = Blueprint("Vizpymzml", __name__, url_prefix='/VizML')

@main.route('/')
def index():
	return render_templat('')

