from flask import Blueprint, request, render_template, \
                  flash, g, session, redirect, url_for
from flask_restful import Api

Vizpymzml = Blueprint("Vizpymzml", __name__, url_prefix='/VizML')
api = Api(Vizpymzml)


@Vizpymzml.route('/')
def index():
	return render_template('VizpyML/index.html')

